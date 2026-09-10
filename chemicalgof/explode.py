from rdkit import Chem
import pandas as pd
from typing import Union, Literal

from .exceptions import InvalidChirality, ConnectorIndex, ExceededBonds
from .gof import DiGraphFrags, FragNode
from .utils import FindProperStereoCenters

import itertools

# return RDKit mol. If invalid reduced graph, it raises error
class Assembler:

    def __init__(self, DiG:DiGraphFrags, strict_chirality:bool = True, verbose = False):
        self.strict_chirality = strict_chirality
        self.verbose = verbose

        self.DiG = DiG

        # adding int on node/fragment atom index to obtain mol atom index
        self.node_add_atom_idx : dict[FragNode:int] = {}

        # atom index and cip label paired to relative node/fragment
        self.node_stereo_atoms : dict[FragNode, tuple[int,str]] = {}

        # atom index and rleative cip label within assembled mol
        self.stereo_atoms : dict[int, int] = {}

        # investigated atoms without a possible chirality
        self.unavailable_chiral_atoms : set[int] = set()

        # atom indexes missing expected stereo cip label
        self.remaining_stereo_atoms : set[int] = set()

        self.single_atom_chiral :set[FragNode] = set()
        self.single_stereo_center : set[FragNode] = set()
        self.psuedo_chiral : set[FragNode] = set()

        self.unsuffixed_chiral :set[FragNode] = set()
        self.other_chiral :set[FragNode] = set()

        self.psuedo_chiral_atoms : dict[FragNode, tuple[int]] = {}

        self.mol : Chem.Mol = None

        self.rank_atoms : dict[int, int] = {}

    def assemble_mol(self):

        UnG = self.DiG.to_undirected()
        node_atom_dummy_label : dict[FragNode, list[tuple[int,int]] ] = {}

        prev_add_idx = 0
        for node in UnG.nodes:

            fragment = node.fragment
            num_atoms = fragment.mol.GetNumAtoms()

            node_atom_dummy_label[node] = []
            self.node_stereo_atoms[node] = list( fragment.chirality.items() )

            self.node_add_atom_idx[node] = prev_add_idx
            prev_add_idx += num_atoms

        for bond_idx, edge in enumerate(UnG.edges, start=1):

            for src, tgt in ( edge, reversed(edge) ):
                src:FragNode
                edge_data = self.DiG.get_edge_data(src, tgt)

                connector_idx = edge_data['aB']

                # TODO
                if not connector_idx in src.fragment.connector_idxs:
                    raise KeyError('Connector index not available')

                node_atom_dummy_label[src].append( (connector_idx, bond_idx) )

                cip_label = edge_data['stereo']
                if cip_label is not None:
                    self.node_stereo_atoms[src].append((connector_idx, cip_label))

        multi_smiles = []
        multi_mols = []

        for node, edges_data in node_atom_dummy_label.items():
            editable_mol = node.fragment.mol
            editable_mol = Chem.EditableMol(editable_mol)

            connections_count : dict[int, int] = {}

            editable_mol.BeginBatchEdit()

            for connector_idx, bond_idx in edges_data:
                new_count = connections_count.get(connector_idx, 0) + 1
                max_count = node.fragment.connector_idxs.get(connector_idx)

                if new_count > max_count:
                    raise ExceededBonds(node.fragment.fragsmiles, connector_idx)

                connections_count[connector_idx] = new_count

                dummy_idx = editable_mol.AddAtom(Chem.AtomFromSmiles(f'[*:{bond_idx}]'))
                editable_mol.AddBond(connector_idx, dummy_idx, Chem.BondType.SINGLE)

            editable_mol.CommitBatchEdit()

            modified_mol = editable_mol.GetMol() # Retrieve the normal molecule object

            for connector_idx, count in connections_count.items():
                atom = modified_mol.GetAtomWithIdx(connector_idx)
                atom.UpdatePropertyCache(strict=False)
                num_exp_HS = atom.GetNumExplicitHs()

                if num_exp_HS > 0:
                    atom.SetNumExplicitHs(num_exp_HS - count) 
                    atom.UpdatePropertyCache()

            # XXX we need mol objects ?
            multi_mols.append(modified_mol)

            # XXX canonical=False because of extended indexes from dummy atoms
            modified_smiles = Chem.MolToSmiles(modified_mol, canonical=False)
            multi_smiles.append(modified_smiles)

        combined_mol = Chem.molzip(Chem.MolFromSmiles('.'.join(multi_smiles)))
        self.mol= combined_mol

        self.rank_atoms.update( dict(
            enumerate(Chem.CanonicalRankAtoms(combined_mol, breakTies=False, includeChirality=False))
        ) )

    def detect_chirality_types(self):

        normal_potential_stereo = set([
            atom_idx for atom_idx, unassigned_str in
            Chem.FindMolChiralCenters(self.mol, force=True, includeUnassigned=True, useLegacyImplementation=True)
        ])

        # XXX useLegacyImplementation=False is necessary
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

        # print(all_potential_stereo, normal_potential_stereo, dependent_potential_stereo)

        for node, stereo_atoms in self.node_stereo_atoms.items():
            if not stereo_atoms:
                continue

            connecting_idxs = set()
            group_rank_atoms = {}

            additive_node_idx = self.node_add_atom_idx[node]
            for connecting_idx, cip_label in stereo_atoms:
                atom_idx = connecting_idx + additive_node_idx

                preset_cip_label = self.stereo_atoms.get(atom_idx)
                if preset_cip_label is not None and cip_label != cip_label:
                    raise KeyError('Double stereo label set on atom index', atom_idx)
                self.stereo_atoms[atom_idx] = cip_label

                # connecting_idxs.append(connecting_idx)
                connecting_idxs.add(atom_idx)

                atom_rank = self.rank_atoms[atom_idx]
                if atom_rank not in group_rank_atoms:
                    group_rank_atoms[atom_rank] = set()

                group_rank_atoms[atom_rank].add(atom_idx)

                if atom_idx not in all_potential_stereo:
                    self.unavailable_chiral_atoms.add(atom_idx)

            fragment = node.fragment

            probab_pseudo = connecting_idxs.intersection(dependent_potential_stereo)
            if probab_pseudo:
                self.psuedo_chiral_atoms[node] = tuple(probab_pseudo)
            elif fragment.num_connector == 1:
                self.single_atom_chiral.add(node)
                continue
            elif len(connecting_idxs) == 1:
                self.single_stereo_center.add(node)
                continue

            # if fragment.group_symmetric_connector:
            #     atom_idxs_fragment = tuple(range( additive_node_idx, fragment.mol.GetNumAtoms() + additive_node_idx ))

            #     equal_ranked_atoms = tuple([ tuple(group) for group in group_rank_atoms.values() if len(group) > 1])

            #     chained_ranked_atoms = tuple(itertools.chain.from_iterable(equal_ranked_atoms))

            #     probab_pseudo = connecting_idxs.difference(chained_ranked_atoms)
            #     if self.verbose:
            #         print(f'{equal_ranked_atoms=}', f'{atom_idxs_fragment=}', f'{group_rank_atoms=}')
            #         print(f'{chained_ranked_atoms=}', f'{connecting_idxs=}', f'{probab_pseudo=}', )

            # if probab_pseudo:
            #     self.psuedo_chiral_atoms.append(tuple(probab_pseudo))


            #     has_pseudo = False

            #     for group in fragment.group_symmetric_connector:
            #         intersection = set(group).intersection(connecting_idxs)
            #         if len(intersection) == 0:
            #             continue
            #         elif len(intersection) == 1:
            #             # print('Pseudo chirality not set correctly.')
            #             # print('Stereo idx:',intersection, 'among', group)
            #             ...

            #         self.psuedo_chiral_atoms.append(
            #             tuple([connecting_idx + additive_node_idx for connecting_idx in intersection])
            #         )

            #         has_pseudo = True

            #     if has_pseudo:
            #         self.psuedo_chiral.add(node)
            #         continue

            elif not fragment.suffix:
                self.unsuffixed_chiral.add(node)
            else:
                self.other_chiral.add(node)

    def initialize_simple_stereocenters(self):

        psuedo_idxs = tuple(itertools.chain.from_iterable(self.psuedo_chiral_atoms.values()))

        # Initializing stereochemistry tag
        for atom_idx in self.stereo_atoms.keys():

            # TODO check unavailable
            if atom_idx in psuedo_idxs or atom_idx in self.unavailable_chiral_atoms:
                continue

            chirA=self.mol.GetAtomWithIdx(atom_idx)
            # chirA.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            chirA.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            self.remaining_stereo_atoms.add(atom_idx)

    def invert_simple_stereocenters_bkp(self):

        atom_idxs = tuple(self.remaining_stereo_atoms)
        # Chem.AssignCIPLabels(self.mol, atomsToLabel=atom_idxs)
        Chem.AssignStereochemistry(self.mol, cleanIt=True, force=True)

        for atom_idx in atom_idxs:

            # TODO
            if self.verbose and atom_idx in self.unavailable_chiral_atoms:
                print('Atom', atom_idx, 'cannot be a stereocenter')

            atom = self.mol.GetAtomWithIdx(atom_idx)

            has_cip = atom.HasProp("_CIPCode") == 1
            if not has_cip:
                if self.verbose:
                    print('Atom', atom_idx, 'has not CIP label')
                continue

            current_cip_label = atom.GetProp("_CIPCode")
            expected_cip_label = self.stereo_atoms[atom_idx]

            if self.verbose and current_cip_label.upper() == expected_cip_label.upper() and current_cip_label != expected_cip_label:
                print('Atom', atom_idx, 'has a incorrect pseudo chirality but correct CIP label')
                print('Expected:',expected_cip_label, 'Current:',current_cip_label)
            elif current_cip_label != expected_cip_label:
                atom.InvertChirality()

            self.remaining_stereo_atoms.discard(atom_idx)

        # Chem.AssignCIPLabels(self.mol, atomsToLabel=atom_idxs)
        Chem.AssignStereochemistry(self.mol, cleanIt=True, force=True)

    def invert_simple_stereocenters(self):

        atom_idxs = tuple(self.remaining_stereo_atoms)
        # Chem.AssignCIPLabels(self.mol, atomsToLabel=atom_idxs)
        Chem.AssignStereochemistry(self.mol, cleanIt=True, force=True)

        current_stereo_centers = FindProperStereoCenters(self.mol, warning=False)
        if self.verbose:
            print(f'{current_stereo_centers=}')

        for atom_idx in atom_idxs:

            # TODO
            if self.verbose and atom_idx in self.unavailable_chiral_atoms:
                print('Atom', atom_idx, 'cannot be a stereocenter')

            current_cip_label = current_stereo_centers.get(atom_idx)
            if current_cip_label is None:
                if self.verbose:
                    print('Atom', atom_idx, 'has not CIP label')
                continue

            expected_cip_label = self.stereo_atoms[atom_idx]

            if self.verbose and current_cip_label.upper() == expected_cip_label.upper() and current_cip_label != expected_cip_label:
                print('Atom', atom_idx, 'has a incorrect pseudo chirality but correct CIP label')
                print('Expected:',expected_cip_label, 'Current:',current_cip_label)
            elif current_cip_label != expected_cip_label:
                atom = self.mol.GetAtomWithIdx(atom_idx)
                atom.InvertChirality()

            self.remaining_stereo_atoms.discard(atom_idx)

        # Chem.AssignCIPLabels(self.mol, atomsToLabel=atom_idxs)
        Chem.AssignStereochemistry(self.mol, cleanIt=True, force=True)

    def force_pseudo_chiral(self):
        chiral_tags_pair = (
                Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        )

        for node, group_pseudo_idxs in self.psuedo_chiral_atoms.items():

            atoms:list[Chem.Atom] = [self.mol.GetAtomWithIdx(atom_idx) for atom_idx in group_pseudo_idxs]
            expected_cip_labels = {atom_idx : self.stereo_atoms[atom_idx] for atom_idx in group_pseudo_idxs}

            additive_node_idx = self.node_add_atom_idx[node]
            fragment_atom_idxs = tuple(range( additive_node_idx, node.fragment.mol.GetNumAtoms() + additive_node_idx ))

            if self.verbose:
                print(f'{group_pseudo_idxs=}', f'{expected_cip_labels=}', f'{fragment_atom_idxs=}')

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

                    # if atom.HasProp("_CIPCode") == 1:
                    #     new_cip_label = atom.GetProp("_CIPCode")
                    # else:
                    #     new_cip_label = '?'
                    new_cip_label = current_stereo_centers.get(atom_idx, '?')

                    new_cip_labels[atom_idx] = new_cip_label
                    matched = new_cip_label.upper() == self.stereo_atoms[atom_idx].upper()

                    solved_pseudo_group.append(matched)
                    # if matched:
                    #     print('Matching', new_cip_label, self.stereo_atoms[atom_idx], 'for atom', atom_idx)

                if self.verbose:
                    print('Current:',new_cip_labels,'Chiral tags:',chiral_names)

                if all(solved_pseudo_group):
                    solved=True
                    break

            if not solved:
                for atom in atoms:
                    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

                # Chem.AssignStereochemistry(self.mol, force=True)
                # Chem.AssignCIPLabels(self.mol, atomsToLabel=group_pseudo_idxs)

                if self.verbose:
                    print('Not solved for', f'{group_pseudo_idxs=}')
                

# XXX return_assembler=True for dev test
def GoF2Mol(DiG:DiGraphFrags, strict_chirality:bool=True, return_assembler=False) -> 'Chem.Mol':
    """Explode reduced graph into rdkit mol object.

    Args:
        DiG (DiGraphFrags): Reduced graph
        strict_chirality (bool, optional): If consider invalid chirality labels provided for atoms. Defaults to True.

    Raises:
        InvalidChirality: If scrict_chirality is True and invalid chiality labels are provided.

    Returns:
        Chem.Mol: RDKit mol object obtained by exploding reduced graph.
    """

    assembler = Assembler(DiG, strict_chirality, verbose=return_assembler)
    assembler.assemble_mol()

    assembler.detect_chirality_types()

    if not assembler.stereo_atoms:
        if return_assembler:
            return assembler
        return assembler.mol
    elif strict_chirality and assembler.unavailable_chiral_atoms:
        raise KeyError('Invalid stereocenters provided:', *assembler.unavailable_chiral_atoms)

    assembler.initialize_simple_stereocenters()
    assembler.invert_simple_stereocenters()

    if assembler.psuedo_chiral_atoms:
        assembler.force_pseudo_chiral()

    Chem.AssignStereochemistry(assembler.mol, cleanIt=True, force=True)
    # Chem.AssignCIPLabels(assembler.mol, atomsToLabel=assembler.stereo_atoms.keys())

    if strict_chirality:
        stereo_centers = FindProperStereoCenters(assembler.mol, warning=False)

        for atom_idx, cip_label in assembler.stereo_atoms.items():
            current_label = stereo_centers.get(atom_idx)
            if current_label is None or current_label.upper() != cip_label.upper():
                chirAtom = assembler.mol.GetAtomWithIdx(atom_idx)
                raise InvalidChirality(chirAtom.GetSymbol(), atom_idx, cip_label)

    if return_assembler:
        return assembler

    return assembler.mol
