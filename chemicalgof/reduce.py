from rdkit import Chem
from .utils import CanonizeFragWithDummies, ClearSmiles, FindProperStereoCenters
from .gof import DiGraphFrags, FragNode, Fragment

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import warnings

# FIXME
# from rdkit.Chem import rdCIPLabeler

class Decompositer:
    # default cleavage pattern. exocyclic single bonds but not beetween charged atoms 
    SINGLEXOCYCLICPATT = '[!$([+1,-1]~[-1,+1])]-&!@[*]'

    def __init__(self, cleavage_smarts:str = SINGLEXOCYCLICPATT):
        self.cleavage_smarts = cleavage_smarts
        self.cleavage_pattern = Chem.MolFromSmarts(cleavage_smarts)

        ## initialized lists (pointed mutables) to be filled ##

        self.fragsMap:list[tuple[int]] = []
        self.fragsIdxs :list[int] = []

        # dicts for mapping atom idxs map[idx_frag][atom_idx] -> atom_idx
        self.mapsFrag2Mol:list[dict[int,int]]=[]
        self.mapsMol2Frag:list[dict[int,int]]=[]

    def fragment(self, mol, bondMatches):
        dumLabels=[(0,0) for _ in bondMatches]
        bonds=[mol.GetBondBetweenAtoms(*atoms).GetIdx() for atoms in bondMatches]

        frags:tuple[str] = Chem.FragmentOnBonds(
            mol,
            addDummies=True,
            bondIndices=bonds, 
            dummyLabels=dumLabels,
        )

        fragsMol:tuple[Chem.Mol] = Chem.GetMolFrags(
            frags,
            asMols=True,
            frags=self.fragsIdxs,
            fragsMolAtomMapping=self.fragsMap
        )

        return fragsMol
    
    def frag_bonds(self, fragsMol, bondMatches : tuple[tuple[int, int]]) -> tuple[tuple[tuple[int,int]]]:
        allBondNeighsFrags = [ [] for _ in range(len(fragsMol))]
        for a,b in bondMatches:
            allBondNeighsFrags[self.fragsIdxs[a]].append((a,b))
            allBondNeighsFrags[self.fragsIdxs[b]].append((b,a))

        return allBondNeighsFrags
    
    def ultimate_smiles(self, fragsMol:tuple[Chem.Mol]):
        # initialize list of cleared (no dummy atoms *) canonized smiles for each fragment
        pureSmis:list[str] = []

        for fMol, fMap in zip(fragsMol, self.fragsMap) :

            Chem.RemoveStereochemistry(fMol)
            # clear dummy atoms, canonize mol fragment to map old idxs with new idxs
            fMol, order = CanonizeFragWithDummies(fMol)

            # mapping idxs
            mapFrag2Mol={v:fMap[k] for k,v in order.items()}
            self.mapsFrag2Mol.append(mapFrag2Mol)

            mapMol2Frag=dict( zip(mapFrag2Mol.values(), mapFrag2Mol.keys()) )
            self.mapsMol2Frag.append(mapMol2Frag)

            # assert to have cleared smiles from mol
            s=ClearSmiles(Chem.MolToSmiles(fMol))

            pureSmis.append(s)

        return pureSmis
    
    def setup_nodes_attributes(self, frag_smiles, allChiralAtoms:dict[int,str], allAtomsInter:list[int]):

        nodes_attributes:list[dict[int,str]] = [] # only chirality attributes. str is R or S linked to atom idx

        for s,mapMol2Frag, mapFrag2Mol in zip(frag_smiles, self.mapsMol2Frag, self.mapsFrag2Mol) :
            # TODO I still dont like to put this here but it's mandatory.
            single_connecting_atom = sum([
                atom.GetTotalNumHs()>0
                for atom in Chem.MolFromSmiles(s).GetAtoms()
            ]) == 1

            # initialize chirality dictionary
            node_attributes = {}
            for _,atom_idx in sorted(mapFrag2Mol.items(), key=lambda x: x[0] ):
                # if frag has only one atom for binding, suffix not include atom idx FIXME in gof traverser when we need to write fragSMILES
                if atom_idx in allChiralAtoms and (atom_idx not in allAtomsInter or single_connecting_atom):
                    node_attributes[ mapMol2Frag[atom_idx] ] = allChiralAtoms[atom_idx]

            nodes_attributes.append(node_attributes)

        return nodes_attributes

def Reduce2GoF(
        # smiles:str = None, mol:Chem.Mol = None, # smiles or mol
        smiles_or_mol :str | Chem.Mol = None,
        /,
        capitalize_legacy:bool = False,
        cleavage_smarts:str = Decompositer.SINGLEXOCYCLICPATT,
        **kwargs, # NOTE employed only for deprecation control
    ) -> DiGraphFrags :
    """Reduce atom-based molecular graph (from smiles or mol object rdkit) into reduced graph fragment-based.

    Args:
        smiles_or_mol (str | Chem.Mol): input as smiles or mol object RDKit.
        capitalize_legacy (bool, optional): if True, pseudo chirality labels (r or s) are forced to be capitalized (as actual chirality labels R or S) and stored as data graph.
        cleavage_smarts (str, optional): SMARTS pattern to employ for reduction (fragmentation) rule. Defaults to exocyclic single bonds fragmentation.

    Raises:
        ValueError: if input is invalid.

    Returns:
        DiGraphFrags: Reduced fragment-based molecular graph.
    """

    if smiles_or_mol is None and not kwargs:
        raise RuntimeError('First positional argument is required: SMILES string or Mol object are allowed.')

    elif smiles_or_mol is None and kwargs:

        valid_old_keys = tuple(set(kwargs.keys()).intersection( ('smiles', 'mol') ))

        if not valid_old_keys:
            KeyError(f'Found invalid kwargs: {", ".join(valid_old_keys)}. Employ only first positional argument.')

        first_valid_key, *_ = valid_old_keys

        warnings.warn(
            "'smiles' and 'mol' are keywords deprecated; employ only first positional argument instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        smiles_or_mol = kwargs[first_valid_key]
    
    ## Providing canonical smiles and then canonical molecule representation
    if  isinstance(smiles_or_mol, Chem.Mol):
        smiles_or_mol = Chem.MolToSmiles(smiles_or_mol)
    elif not isinstance(smiles_or_mol, str):
        raise TypeError('Input Error : SMILES string or Chem.Mol object is required.')

    smiles_or_mol = Chem.CanonSmiles(smiles_or_mol)
    mol = Chem.MolFromSmiles(smiles_or_mol)
    
    obj = Decompositer(cleavage_smarts)

    bondMatches:tuple[tuple[int,int]] = mol.GetSubstructMatches( obj.cleavage_pattern )

    # optical stereochemical data
    allChiralAtoms = FindProperStereoCenters(mol)

    if capitalize_legacy and allChiralAtoms:
        allChiralAtoms = {atom_idx: cip_label.upper() for atom_idx,cip_label in allChiralAtoms.items()}

    if not bondMatches: # 1 single fragment -> no edges within graph.
        fragsMol = [mol]
        frag_smiles = [Chem.MolToSmiles(mol)]
        frag_bonds = []
        nodes_attributes = [ {a:label for a,label in sorted(allChiralAtoms.items(), key=lambda x: x[0] )} ]
    else:
        fragsMol = obj.fragment(mol, bondMatches)

        # outer list follows frag idxs
        # inner lists : [ idxAtom_1frag, idxAtom_2frag  ], [ idxAtom_2frag, idxAtom_1frag  ]
        frag_bonds = obj.frag_bonds(fragsMol, bondMatches)

        frag_smiles = obj.ultimate_smiles(fragsMol)

        # store each atoms involved in bonds. It needs to label only chiral atom within fragment, not for edges
        from itertools import chain
        allAtomsInter = set(chain(*bondMatches))

        nodes_attributes = obj.setup_nodes_attributes(frag_smiles, allChiralAtoms, allAtomsInter)

    # initialize directed reduced graph

    #initialize nodes
    fragments = tuple([
        Fragment(smiles=node_smiles, chirality=chirality) 
        for node_smiles,chirality in zip(frag_smiles, nodes_attributes)
    ])

    nodes = tuple([FragNode(fragment) for fragment in fragments])

    # add nodes and their attributes (loaded from node class attributes)
    diG=DiGraphFrags()
    diG.add_nodes_from(nodes)

    # add edges, if there are ...
    for    mapMol2Frag,  node,        fNeig,        in \
    zip(   obj.mapsMol2Frag, nodes, frag_bonds):

        # setup edges in directed graph !
        for master_atom_idx,neighbor_frag_idx in fNeig:
            neigh_node=nodes[obj.fragsIdxs[neighbor_frag_idx]]
            connecting_atom_idx=mapMol2Frag[master_atom_idx]
            stereo = allChiralAtoms.get(master_atom_idx) if node.fragment.num_connector > 1 else None
            diG.add_edge(node, neigh_node, aB=connecting_atom_idx, stereo=stereo )

    return diG

