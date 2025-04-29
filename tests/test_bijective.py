import pytest
from chemicalgof import encode, decode
from rdkit import Chem


def _fast_check(smiles: str) -> str:
    fragsmiles = encode(smiles)
    decoded = decode(fragsmiles, strict_chirality=True)
    return Chem.CanonSmiles(decoded)


@pytest.mark.parametrize(
    "smiles",
    [
        # Capital chirality dependency
        'O=C(O[C@H]1C[C@H]2CC[C@@H](C1)N2C[C@H](O)c1ccccc1)c1ccccc1',

        # Pseudo chirality dependency
        'CCCCc1cn([C@H]2[C@H](C)CCC[C@@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1',

        # Double pseudo chirality
        'Cc1cc2c(cc1Cc1ccc(C(=O)NC[C@H]3CC[C@H](C(N)=O)CC3)o1)C(C)(C)CCC2(C)C'
    ],
    ids=[
        "capital_chirality_dependency",
        "pseudo_chirality_dependency",
        "double_pseudo_chirality"
    ]
)
def test_bijective_smiles(smiles):
    decoded = _fast_check(smiles)
    assert smiles == decoded
