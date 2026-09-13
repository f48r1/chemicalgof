from rdkit import Chem

from .exceptions import InvalidChirality, ConnectorIndex, ExceededBonds, CollidedMultiChirality
from .gof import DiGraphFrags, FragNode, Fragment
from .utils import FindProperStereoCenters

import itertools

# return RDKit mol. If invalid reduced graph, it raises error
class Assembler:

    def __init__(self, DiG:DiGraphFrags, strict_chirality:bool = True):
        self.strict_chirality = strict_chirality

        self.DiG = DiG

        # adding int on node/fragment atom index to obtain mol atom index
        self.node_add_atom_idx : dict[FragNode:int] = {}

        # atom index and cip label paired to relative node/fragment
        self.node_stereo_atoms : dict[FragNode, tuple[int,str]] = {}

        # atom index and rleative cip label within assembled mol
        self.stereo_atoms : dict[int, int] = {}

        # investigated atoms without a possible chirality
        self.unavailable_chiral_atoms : set[int] = set()

        self.psuedo_chiral_atoms : dict[FragNode, tuple[int]] = {}

        self.mol : Chem.Mol = None

    def assemble_mol(self):

        UnG = self.DiG.to_undirected()

        editable_mol = Chem.Mol()
        num_expl_Hs_connectors = {}

        prev_add_idx = 0
        for node in UnG.nodes:

            fragment:Fragment = node.fragment
            frag_mol = fragment.mol
            num_atoms = frag_mol.GetNumAtoms()

            self.node_stereo_atoms[node] = []

            self.node_add_atom_idx[node] = prev_add_idx

            for atom_idx in fragment.connector_idxs.keys():
                atom = frag_mol.GetAtomWithIdx(atom_idx)

                num_exp_HS = atom.GetNumExplicitHs()
                if num_exp_HS > 0:
                    num_expl_Hs_connectors[atom_idx + prev_add_idx] = num_exp_HS

            editable_mol = Chem.CombineMols(editable_mol, frag_mol)

            prev_add_idx += num_atoms

        connections_count : dict[int, int] = {}

        editable_mol = Chem.EditableMol(editable_mol)
        editable_mol.BeginBatchEdit()

        for edge in UnG.edges:

            edge:tuple[FragNode, FragNode]

            atom_idxs_src2tgt = []

            for src, tgt in ( edge, reversed(edge) ):
                edge_data = self.DiG.get_edge_data(src, tgt)

                connector_idx_frag = edge_data['aB']
                fragment = src.fragment

                max_count = fragment.connector_idxs.get(connector_idx_frag)
                if max_count is None:
                    raise ConnectorIndex(fragment.fragsmiles, connector_idx_frag)
                
                connector_idx_mol = connector_idx_frag + self.node_add_atom_idx[src]
                current_count = connections_count.get(connector_idx_mol, 0) + 1 

                if current_count > max_count:
                    raise ExceededBonds(fragment.fragsmiles, connector_idx_frag)

                connections_count[connector_idx_mol] = current_count

                cip_label = edge_data.get('stereo')
                if cip_label is not None:
                    self.node_stereo_atoms[src].append((connector_idx_frag, cip_label))

                atom_idxs_src2tgt.append(connector_idx_mol)

            editable_mol.AddBond(*atom_idxs_src2tgt, Chem.BondType.SINGLE)

        editable_mol.CommitBatchEdit()

        editable_mol = editable_mol.GetMol() # Retrieve the normal molecule object

        for connector_idx, count in connections_count.items():

            num_exp_HS = num_expl_Hs_connectors.get(connector_idx, 0)
            if not num_exp_HS:
                continue
            elif count >= num_exp_HS:
                obtained_num_exp_Hs = 0
            else:
                obtained_num_exp_Hs = num_exp_HS - count

            atom = editable_mol.GetAtomWithIdx(connector_idx)
            atom.SetNumExplicitHs(obtained_num_exp_Hs) 

            atom.UpdatePropertyCache()

        self.mol = editable_mol

    def detect_chirality_types(self):

        normal_potential_stereo = set([
            atom_idx for atom_idx, unassigned_str in
            Chem.FindMolChiralCenters(self.mol, force=True, includeUnassigned=True, useLegacyImplementation=True)
        ])

        try:
            all_potential_stereo = set([
                atom_idx for atom_idx, unassigned_str in
                Chem.FindMolChiralCenters(self.mol, force=True, includeUnassigned=True, useLegacyImplementation=False)
            ])

        except:
            all_potential_stereo = set()

        if not all_potential_stereo:
            all_potential_stereo = normal_potential_stereo.copy()

        dependent_potential_stereo = all_potential_stereo.difference(normal_potential_stereo)

        for node, stereo_atoms in self.node_stereo_atoms.items():
            fragment = node.fragment

            if not stereo_atoms and not fragment.chirality:
                continue

            additive_node_idx = self.node_add_atom_idx[node]
            all_connecting_idx_mol = set()

            def check_connecting_stereo(iterable_atom_label:list[tuple[int,str]], ignore=False):

                for connecting_idx_frag, cip_label in iterable_atom_label:
                    connecting_idx_mol = connecting_idx_frag + additive_node_idx

                    preset_cip_label = self.stereo_atoms.get(connecting_idx_mol)

                    if connecting_idx_mol in self.unavailable_chiral_atoms:
                        continue

                    if preset_cip_label is not None and cip_label.upper() != preset_cip_label.upper():
                        if self.strict_chirality:
                            raise CollidedMultiChirality(fragment.mol.GetAtomWithIdx(connecting_idx_frag).GetSymbol(), connecting_idx_mol)
                        elif not ignore:
                            self.stereo_atoms.pop(connecting_idx_mol)
                            all_connecting_idx_mol.discard(connecting_idx_mol)
                            self.unavailable_chiral_atoms.add(connecting_idx_mol)

                    elif connecting_idx_mol not in all_potential_stereo:
                        if self.strict_chirality:
                            raise InvalidChirality(fragment.mol.GetAtomWithIdx(connecting_idx_frag).GetSymbol(), connecting_idx_mol, cip_label)
                        
                        self.unavailable_chiral_atoms.add(connecting_idx_mol)
                    elif preset_cip_label is None:
                        all_connecting_idx_mol.add(connecting_idx_mol)
                        self.stereo_atoms[connecting_idx_mol] = cip_label

            check_connecting_stereo(stereo_atoms)
            check_connecting_stereo(fragment.chirality.items(), ignore=True)

            probab_pseudo = all_connecting_idx_mol.intersection(dependent_potential_stereo)
            if probab_pseudo:
                self.psuedo_chiral_atoms[node] = tuple(probab_pseudo)

    def label_independent_stereocenters(self):

        psuedo_idxs = tuple(itertools.chain.from_iterable(self.psuedo_chiral_atoms.values()))

        # Initializing stereochemistry tag
        for atom_idx in self.stereo_atoms.keys():

            if atom_idx in psuedo_idxs:
                continue

            atom=self.mol.GetAtomWithIdx(atom_idx)

            # atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)

        # Chem.AssignCIPLabels(self.mol, atomsToLabel=atom_idxs)
        Chem.AssignStereochemistry(self.mol, cleanIt=True, force=True)

        current_stereo_centers = FindProperStereoCenters(self.mol, warning=False)
        for atom_idx in self.stereo_atoms.keys():

            if atom_idx in psuedo_idxs:
                continue

            current_cip_label = current_stereo_centers.get(atom_idx)
            if current_cip_label is None:
                if self.strict_chirality:
                    raise InvalidChirality()
                
                self.unavailable_chiral_atoms.add(atom_idx)
                continue

            expected_cip_label = self.stereo_atoms[atom_idx]

            if current_cip_label.upper() != expected_cip_label.upper() :
                atom = self.mol.GetAtomWithIdx(atom_idx)
                atom.InvertChirality()

        # Chem.AssignCIPLabels(self.mol, atomsToLabel=atom_idxs)
        Chem.AssignStereochemistry(self.mol, cleanIt=True, force=True)

    def force_dependent_stereocenters(self):
        chiral_tags_pair = (
                Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        )

        for node, group_pseudo_idxs in self.psuedo_chiral_atoms.items():

            atoms:list[Chem.Atom] = [self.mol.GetAtomWithIdx(atom_idx) for atom_idx in group_pseudo_idxs]

            solved = False
            for chiral_tags in itertools.product(chiral_tags_pair, repeat=len(atoms)):

                chiral_names = []

                for element_idx, chiral_tag in enumerate(chiral_tags):
                    atom_idx = group_pseudo_idxs[element_idx]
                    atom = atoms[element_idx]

                    atom.SetChiralTag(chiral_tag)
                    chiral_names.append(str(chiral_tag))

                Chem.AssignStereochemistry(self.mol, cleanIt=False, force=True)
                # Chem.AssignCIPLabels(self.mol, atomsToLabel=self.stereo_atoms.keys())
                current_stereo_centers = FindProperStereoCenters(self.mol, warning=False)

                solved_pseudo_group = []
                new_cip_labels = {}

                for element_idx in range(len(group_pseudo_idxs)):
                    atom_idx = group_pseudo_idxs[element_idx]
                    atom = atoms[element_idx]

                    new_cip_label = current_stereo_centers.get(atom_idx, '?')

                    new_cip_labels[atom_idx] = new_cip_label
                    matched = new_cip_label.upper() == self.stereo_atoms[atom_idx].upper()

                    solved_pseudo_group.append(matched)

                if all(solved_pseudo_group):
                    solved=True
                    break

            if not solved:
                if self.strict_chirality:
                    raise InvalidChirality()
                
                for atom in atoms:
                    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

                # Chem.AssignStereochemistry(self.mol, force=True)
                # Chem.AssignCIPLabels(self.mol, atomsToLabel=group_pseudo_idxs)
                
def GoF2Mol(DiG:DiGraphFrags, strict_chirality:bool=True) -> 'Chem.Mol':
    """Explode reduced graph into rdkit mol object.

    Args:
        DiG (DiGraphFrags): Reduced graph
        strict_chirality (bool, optional): If consider invalid chirality labels provided for atoms. Defaults to True.

    Raises:
        InvalidChirality: If scrict_chirality is True and invalid chiality labels are provided.

    Returns:
        Chem.Mol: RDKit mol object obtained by exploding reduced graph.
    """

    assembler = Assembler(DiG, strict_chirality)
    assembler.assemble_mol()

    assembler.detect_chirality_types()

    if not assembler.stereo_atoms:
        return assembler.mol

    assembler.label_independent_stereocenters()

    if assembler.psuedo_chiral_atoms:
        assembler.force_dependent_stereocenters()

        Chem.AssignStereochemistry(assembler.mol, cleanIt=True, force=True)
        # Chem.AssignCIPLabels(assembler.mol, atomsToLabel=assembler.stereo_atoms.keys())

    if strict_chirality:
        stereo_centers = FindProperStereoCenters(assembler.mol, warning=False)

        for atom_idx, cip_label in assembler.stereo_atoms.items():
            current_label = stereo_centers.get(atom_idx)
            if current_label is None or current_label.upper() != cip_label.upper():
                chirAtom = assembler.mol.GetAtomWithIdx(atom_idx)
                raise InvalidChirality(chirAtom.GetSymbol(), atom_idx, cip_label)

    return assembler.mol
