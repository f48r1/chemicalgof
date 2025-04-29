import os
import pytest
import pandas as pd
from random import choice
from rdkit import Chem

from chemicalgof import (
    encode,
    decode,
)

# === Fixtures ===

@pytest.fixture(scope="session")
def dataset():
    dataset_path = os.path.join(os.path.dirname(__file__), "../data/test.csv")
    return pd.read_csv(dataset_path, header=None).squeeze()

@pytest.fixture
def random_bug_node(dataset):
    bug_nodes_idxs = [113284, 113285, 113286]
    idx = choice(bug_nodes_idxs)
    return dataset.loc[idx]

@pytest.fixture
def random_bug_chirality(dataset):
    bug_chirality_idxs = [74558]
    idx = choice(bug_chirality_idxs)
    return dataset.loc[idx]

# === Helper ===

def _fast_check(smiles: str) -> str:
    fragsmiles = encode(smiles)
    decoded = decode(fragsmiles, strict_chirality=True)
    return Chem.CanonSmiles(decoded)

# === Test functions ===

def test_nodes_overcount(random_bug_node):
    smiles = random_bug_node

    with pytest.raises(Exception) as exc_info:
        decoded = _fast_check(smiles)

    assert "more than 100000nodes found" in str(exc_info.value)

def test_chirality_unrecognized(random_bug_chirality):
    smiles = random_bug_chirality
    mol = Chem.MolFromSmiles(smiles)

    merged_stereo = dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)) # NOTE bug: we need either legacy=True and False
    merged_stereo = merged_stereo | dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=True))
		
    default_stereo = dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=False))

    keys_check = set(list(merged_stereo.keys())) == set(list(default_stereo.keys()))

    assert not keys_check

    decoded = _fast_check(smiles)
    smiles_canonical = Chem.CanonSmiles(smiles)

    assert decoded != smiles_canonical
