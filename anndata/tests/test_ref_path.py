import pytest
import numpy as np
import pandas as pd
import scipy.sparse as ssp

from anndata import AnnData
from anndata._core.ref_path import RefPath, split_paths, resolve_path


@pytest.fixture
def adata():
    return AnnData(
        X=np.array([[0, 1, 2, 3], [4, 5, 6, 7]]),
        layers=dict(unspliced=np.array([[0.1, 1.2, 2.3, 3.4], [4.5, 5.6, 6.7, 7.8]])),
        obs=dict(obs_names=["Cell1", "Cell2"], group=["batch1", "batch2"]),
        obsm=dict(
            X_pca=np.array([[0.2, 0.4, 0.6], [0.1, 0.3, 0.5]]),
            protein=pd.DataFrame(dict(CD14=[2.1, 3.2])),
        ),
        obsp=dict(neighbors_distances=ssp.csr_matrix([[0, 0.4], [0.3, 0]])),
        var=dict(var_names=[f"Gene{i}" for i in "WXYZ"], mito=[0.4, 0.3, 0.2, 0.1]),
        # varm=dict(...),
        varp=dict(cors=np.array([(5, 4, 3, 2)] * 4)),
        raw=dict(
            X=np.array([[0, 11, 22, 33, 88], [44, 55, 66, 77, 99]]),
            var=dict(
                var_names=[f"Gene{i}" for i in "WXYZA"],
                symbol=[f"Symb{i}" for i in range(1, 6)],
            ),
        ),
    )


paths = [
    (("layers", "unspliced", "obs", "Cell1"), "var", [0.1, 1.2, 2.3, 3.4]),
    (("obs", "group"), "obs", ["batch1", "batch2"]),
    (("obsm", "X_pca", 1), "obs", [0.4, 0.3]),
    (("obsm", "protein", "CD14"), "obs", [2.1, 3.2]),
    (("obsp", "neighbors_distances", "Cell2", 0), "obs", [0.3, 0]),
    (("var", "mito"), "var", [0.4, 0.3, 0.2, 0.1]),
    # (("varm", "", ""), "var", []),
    (("varp", "cors", "GeneY", 1), "var", [3] * 4),
    (("raw", "X", "var", "GeneA"), "obs", [88, 99]),
    # TODO: is var correct here? It’ll return more variables …
    (("raw", "var", "symbol"), "var", ["Symb1", "Symb2", "Symb3", "Symb4", "Symb5"]),
]


@pytest.mark.parametrize("args,dim,expected", paths, ids=repr)
def test_dim(args, dim, expected):
    assert RefPath(*args).dim == dim


@pytest.mark.parametrize("args,dim,expected", paths, ids=repr)
def test_get_vector(args, dim, expected, adata):
    rp = RefPath(*args)
    vec = rp.get_vector(adata)
    assert isinstance(vec, np.ndarray)
    assert len(vec.shape) == 1
    assert vec.tolist() == expected


@pytest.mark.parametrize(
    "rp_code", ["RefPath('obs', 'group')", "RefPath('varp', 'neighbors', 'Cell5', 1)"]
)
def test_repr(rp_code):
    rp = eval(rp_code)
    assert repr(rp) == rp_code


@pytest.mark.parametrize(
    "spec, resolved",
    [
        (("X", "var", "GeneX"), ("layers", "X", "var", "GeneX")),
        ("op/foo/bar/0", ("obsp", "foo", "bar", 0)),
        (("rX", "var", "GeneA"), ("raw", "X", "var", "GeneA")),
    ],
    ids=repr,
)
def test_parse(spec, resolved):
    assert RefPath.parse(spec) == RefPath(*resolved)


@pytest.mark.parametrize(
    "spec,err_regex",
    [
        (("raw",), r"No path specified\."),
        (("obsm", "X_pca", None), r"None in path \('X_pca', None\).*\{int, str\}\."),
        (("layers", "X", 1, "G1"), r"path\[1\]=1.*not one of \{'obs', 'var'\}\."),
        ("X", r"required \(3\) ≠ got \(1\) in path \('X',\) for attr layers\."),
        ("f/XY", r"Unknown attr name or short code 'f'\."),
        ("op/foo/bar/notAnInt", r"Invalid.*obsp path: 'notAnInt'\."),
    ],
    ids=repr,
)
def test_parse_failures(spec, err_regex):
    with pytest.raises(ValueError, match=err_regex):
        RefPath.parse(spec)


# anndata part (resolving paths)


def test_alias():
    pass  # TODO


@pytest.mark.parametrize(
    "multipath,reference",
    [
        ([("obs", "A"), ("X", "var", "G1")], [("obs", "A"), ("X", "var", "G1")]),
        (("obs", ["A", "B"]), [("obs", "A"), ("obs", "B")]),
        (
            (["obs", "var"], ["A", "B"]),
            [("obs", "A"), ("obs", "B"), ("var", "A"), ("var", "B")],
        ),
        (
            ["G1", ("rX", "var", ["G4", "G5"])],
            ["G1", ("rX", "var", "G4"), ("rX", "var", "G5")],
        ),
    ],
    ids=repr,
)
def test_split_paths(multipath, reference):
    got = list(split_paths(multipath))
    assert got == reference


@pytest.mark.parametrize(
    "short_path,resolved",
    [
        # single strings should be found in {obs,var}{.columns,_names}
        # but not in {obs,var}m DataFrame columns (too deep and obscure)
        ("Cell1", ("layers", "X", "obs", "Cell1")),
        ("group", ("obs", "group")),
        # keys with subpaths should be found in layers, {obs,var}{p,m}
        ("X_pca/1", ("obsm", "X_pca", 1)),
        ("unspliced/GeneY", ("layers", "unspliced", "var", "GeneY")),
        # {obs,var}p should default to axis 0
        (("neighbors_distances", "Cell2"), ("obsp", "neighbors_distances", "Cell2", 0)),
    ],
    ids=repr,
)
def test_resolve(adata, short_path, resolved):
    assert resolve_path(adata, short_path) == RefPath(*resolved)
