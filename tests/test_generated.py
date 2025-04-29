import pytest
from chemicalgof import decode
from rdkit import Chem


@pytest.mark.parametrize(
    "sampled",
    [
        # Canonization is required beacause RDKit assign chirality to sp2 catbon atoms
        ['C', '<12R>', 'O=C1NCCCc2cccc(c2)CCCOCc2cccc1c2', '<20R>', '<16>', '(', 'C', ')', '<11R>', '(', '<6>', 'c1ccc2c(c1)CCC21CCNCC1', ')', '<22S>', '(', 'O', ')', 'O', 'C'],

    ],
    ids=[
        "required_canonization",
    ]
)

def test_decoding_sampled(sampled):
    decoded = decode(sampled)
    canonical = Chem.CanonSmiles(decoded)
    assert decoded != canonical