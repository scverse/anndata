import pytest

from anndata._core.ref_path import RefPath


init_args = [
    ("layers", "foo", "obs", "Cell1"),
    ("obs", "group"),
    ("obsm", "X_pca", 1),
    ("obsp", "neighbors_distances", "Cell5", 1),
    ("var", "percent_mito"),
    ("varm", "protein", "CD14"),
    ("varp", "cors", "GeneY", 0),  # no idea what to put here
    ("raw", "X", "var", "GeneX"),
]


@pytest.mark.parametrize("args", init_args)
def test_instantiate(args):
    RefPath(*args)


@pytest.mark.parametrize(
    "rp_code", ["RefPath('obs', 'group')", "RefPath('varp', 'neighbors', 'Cell5', 1)"]
)
def test_repr(rp_code):
    rp = eval(rp_code)
    assert repr(rp) == rp_code


def test_get_vector():
    pass
