import pytest
from chemicalgof import encode

# === Fixtures ===

@pytest.fixture(scope="session")
def dataset():
    import pandas as pd
    from pathlib import Path

    dataset_path = Path(__file__).parent.parent / "data" / "test.csv"
    return pd.read_csv(dataset_path, header=None).squeeze()

@pytest.fixture
def random_warning_node(dataset):
    warning_nodes_idxs = [113284, 113285, 113286]
    subset = dataset.loc[warning_nodes_idxs]
    return subset.sample(n=1, random_state=None).item()

# === Test functions ===

def test_nodes_overcount(random_warning_node):
    smiles = random_warning_node

    with pytest.warns(RuntimeWarning) as warning_info:
        fragsmiles = encode(smiles)

    assert "more than 100000nodes found" in str(warning_info[0].message)
