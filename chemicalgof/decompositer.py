from rdkit import Chem
from .utils import CanonizeFragWithDummies, ClearSmiles, GetPotAtomLinkers
from .classes import DiGraphFrags, FragNode
from rdkit import RDLogger
# from rdkit.Chem import rdCIPLabeler
import numpy as np
RDLogger.DisableLog('rdApp.*')


# default cleavage pattern. exocyclic single bonds but not beetween charged atoms 
SINGLEXOCYCLICPATT = '[!$([+1,-1]~[-1,+1])]-&!@[*]'

def getOneSingleFragSmiles(mol):

    ## TODO
    ## Check if stereochemistry removal leads to different atom indexing
    allChiralAtoms=dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=True ))
    
    Chem.RemoveStereochemistry(mol)
    # clear dummy atoms, canonize mol fragment to map old idxs with new idxs
    mol, order = CanonizeFragWithDummies(mol)

    s = Chem.MolToSmiles(mol)

    # initialize stereochemical suffix
    suff=""
    for a,RS in sorted(allChiralAtoms.items(), key=lambda x: x[0] ):
        # if frag has only one atom for binding, suffix not include atom idx
        if len(GetPotAtomLinkers(s))==1:
            suff+=RS
        # single fragment never includes linker atoms !
        # elif a not in allAtomsInter:
        else:
            suff+=str(a)+RS
    # if stereochemical is provided, suffix is added
    if suff:
        return s+"|"+suff
    
    return s

def Mol2GoF (mol:Chem.Mol, pattBonds:str = SINGLEXOCYCLICPATT, ) -> DiGraphFrags | None:
    ## Providing canonical smiles and then canonical molecule representation
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    bondMatches = mol.GetSubstructMatches( Chem.MolFromSmarts(pattBonds) )
    if not bondMatches:
        s = getOneSingleFragSmiles(mol)

        # define directed reduced graph
        diG=DiGraphFrags()
        diG.add_node(FragNode.fromSmiles(s))
        return diG

    bondMatches = np.array(bondMatches)

    # store each atoms involved in bonds
    allAtomsInter = np.unique( bondMatches.flatten() )
    
    dumLabels=[(0,0) for _ in bondMatches]
    bonds=[mol.GetBondBetweenAtoms(*b).GetIdx() for b in bondMatches.tolist()]  
    frags = Chem.FragmentOnBonds(mol, 
                                    addDummies=True,
                                    bondIndices=bonds, 
                                    dummyLabels=dumLabels,
                                )
    fragsMap=[]
    fragsIdxs=[]
    fragsMol=Chem.GetMolFrags(frags,asMols=True,frags=fragsIdxs,fragsMolAtomMapping=fragsMap)
    
    # optical stereochemical data
    allChiralAtoms=dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=True ))
    

    # initialize list of lists to be filled of neighbours data:
    # outer list follows frag idxs
    # inner lists : [ idxAtom_1frag, idxAtom_2frag  ], [ idxAtom_2frag, idxAtom_1frag  ]
    allBondNeighsFrags=[[] for _ in range(len(fragsMol))]
    for a,b in bondMatches:
        allBondNeighsFrags[fragsIdxs[a]].append((a,b))
        allBondNeighsFrags[fragsIdxs[b]].append((b,a))

    # initialize list of cleared (no dummy atoms *) canonized smiles for each fragment
    pureSmis = []

    # dicts for mapping atom idxs 
    mapsFrag2Mol=[]
    mapsMol2Frag=[]

    for fMol, fMap in zip(fragsMol, fragsMap) :

        Chem.RemoveStereochemistry(fMol)
        # clear dummy atoms, canonize mol fragment to map old idxs with new idxs
        fMol, order = CanonizeFragWithDummies(fMol)

        # mapping idxs
        mapFrag2Mol={v:fMap[k] for k,v in order.items()}
        mapsFrag2Mol.append(mapFrag2Mol)

        mapMol2Frag=dict( zip(mapFrag2Mol.values(), mapFrag2Mol.keys()) )
        mapsMol2Frag.append(mapMol2Frag)

        # assert to have cleared smiles from mol
        s=ClearSmiles(Chem.MolToSmiles(fMol))

        # initialize stereochemical suffix
        suff=""
        for a in [v for _,v in sorted(mapFrag2Mol.items(), key=lambda x: x[0] ) if v in allChiralAtoms]:
            # if frag has only one atom for binding, suffix not include atom idx
            if len(GetPotAtomLinkers(s))==1:
                suff+=allChiralAtoms[a]
            # else it includes atom idxs
            elif a not in allAtomsInter:
                suff+=str(mapMol2Frag[a])+allChiralAtoms[a]
        # if stereochemical is provided, suffix is added
        if suff:
            pureSmis.append ( s+"|"+suff )
        else:
            pureSmis.append ( s )

    ##
    # ad MV, che la strada mi illumina ..
    ##

    # define actual nodes belonging to reduced graph.
    nodes=[FragNode.fromSmiles(x) for x in pureSmis]

    # define directed reduced graph
    diG=DiGraphFrags()
    diG.add_nodes_from(nodes)

    # zipping node info for each fragment
    for    mapMol2Frag,  node,        fNeig,        in \
    zip(   mapsMol2Frag, nodes, allBondNeighsFrags):

        # setup node in directed graph !
        for a,nn in fNeig:
            neigh=nodes[fragsIdxs[nn]]
            aB=mapMol2Frag[a]
            if a in allChiralAtoms:
                diG.add_edge(node, neigh, aB=aB, stereo=allChiralAtoms[a] )
            else:
                diG.add_edge(node, neigh, aB=aB)

    return diG

def Smiles2GoF (smiles:str, pattBonds:str = SINGLEXOCYCLICPATT, ) -> DiGraphFrags | None:
    return Mol2GoF(Chem.MolFromSmiles(smiles), pattBonds)