from .reduce import Reduce2GoF
from .write import GoF2fragSMILES
from .explode import GoF2Mol
from .parse import fragSMILES2GoF, Sequence2GoF, split

def encode(
    smiles:str,
    canonical:bool=True,
    random:bool=False,
    capitalize_chirality=True,
) -> str :
    """Convert SMILES string into fragSMILES representation. Reduced graph is considered as an intermediate of this function but it is not shown.

    Args:
        smiles (str): _description_
        canonical (bool, optional): Traversing intermediate reduced graph by canonical way. Defaults to True.
        random (bool, optional): Reduced intermediate graph is traversed randomly. Defaults to False.
        capitalize_chirality (bool, optional): If consider pseudo-chirality (r or s labels) as a actual chirality (R or S labels). Defaults to True.

    Returns:
        str: fragSMILES representation. Then string can be splitted by function provided by this package.
    """
    DiG = Reduce2GoF(smiles=smiles, capitalize_legacy=capitalize_chirality)
    fragsmiles = GoF2fragSMILES(DiG, canonize=canonical, random=random)

    return fragsmiles
    
def decode(
    fragsmiles:str | list[str],
    strict_chirality:bool = True,
) -> str :
    """Convert fragSMILES string (if string) or fragSMILES tokenized sequence (if list of strings) to relative SMILES.
    
    Args:
        fragsmiles (str | list[str]): string of fragSMILES representation or fragSMILES tokenized sequence (usefull if coversion involves generated sequence by models).
        strict_chirality (bool, optional): If take in account invalid assigned chirality label. Raise error when invalid labels are provided. Defaults to True.

    Returns:
        str: SMILES string. Tip: Canonize it for correct representation and compare invalid chiral atoms.
    """
    from rdkit import Chem
    if type(fragsmiles) is str:
        DiG = fragSMILES2GoF(fragsmiles)
    else:
        DiG = Sequence2GoF(fragsmiles)

    mol = GoF2Mol(DiG, strict_chirality=strict_chirality)
    smiles = Chem.MolToSmiles(mol)
    # smiles = Chem.CanonSmiles(smiles) # [x] Canonization is not preferred because of bug about chirality: it's still expected for aromatic and sp2 carbon atoms. If you canonize returned SMILES, sanification can be done on it!
    return smiles