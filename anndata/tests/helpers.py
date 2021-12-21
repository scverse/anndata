from string import ascii_letters
from typing import Tuple
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from anndata import AnnData


def gen_vstr_recarray(m, n, dtype=None):
    size = m * n
    lengths = np.random.randint(3, 5, size)
    letters = np.array(list(ascii_letters))
    gen_word = lambda l: "".join(np.random.choice(letters, l))
    arr = np.array([gen_word(l) for l in lengths]).reshape(m, n)
    return pd.DataFrame(arr, columns=[gen_word(5) for i in range(n)]).to_records(
        index=False, column_dtypes=dtype
    )


def gen_typed_df(n, index=None):
    # TODO: Think about allowing index to be passed for n
    letters = np.fromiter(iter(ascii_letters), "U1")
    if n > len(letters):
        letters = letters[: n // 2]  # Make sure categories are repeated
    return pd.DataFrame(
        dict(
            cat=pd.Categorical(np.random.choice(letters, n)),
            cat_ordered=pd.Categorical(np.random.choice(letters, n), ordered=True),
            int64=np.random.randint(-50, 50, n),
            float64=np.random.random(n),
            uint8=np.random.randint(255, size=n, dtype="uint8"),
        ),
        index=index,
    )


def gen_typed_df_t2_size(m, n, index=None, columns=None) -> pd.DataFrame:
    s = 0
    df = pd.DataFrame()
    new_vals = gen_typed_df(m)
    while s < (n / new_vals.shape[1]):
        new_vals = gen_typed_df(m, index=index)
        new_vals.columns = new_vals.columns + "_" + str(s)
        df[new_vals.columns] = new_vals
        s += 1
    df = df.iloc[:m, :n].copy()
    if columns is not None:
        df.columns = columns
    return df


# TODO: Use hypothesis for this?
def gen_adata(
    shape: Tuple[int, int],
    X_type=sparse.csr_matrix,
    X_dtype=np.float32,
    # obs_dtypes,
    # var_dtypes,
    obsm_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame),
    varm_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame),
    layers_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame),
) -> AnnData:
    """\
    Helper function to generate a random AnnData for testing purposes.

    Note: For `obsm_types`, `varm_types`, and `layers_types` these currently
    just filter already created objects.
    In future, these should choose which objects are created.

    Params
    ------
    shape
        What shape you want the anndata to be.
    X_type
        What kind of container should `X` be? This will be called on a randomly
        generated 2d array.
    X_dtype
        What should the dtype of the `.X` container be?
    obsm_types
        What kinds of containers should be in `.obsm`?
    varm_types
        What kinds of containers should be in `.varm`?
    layers_types
        What kinds of containers should be in `.layers`?
    """
    M, N = shape
    obs_names = pd.Index(f"cell{i}" for i in range(shape[0]))
    var_names = pd.Index(f"gene{i}" for i in range(shape[1]))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    # For #147
    obs.rename(columns=dict(cat="obs_cat"), inplace=True)
    var.rename(columns=dict(cat="var_cat"), inplace=True)

    if X_type is None:
        X = None
    else:
        X = X_type(np.random.binomial(100, 0.005, (M, N)).astype(X_dtype))
    obsm = dict(
        array=np.random.random((M, 50)),
        sparse=sparse.random(M, 100, format="csr"),
        df=gen_typed_df(M, obs_names),
    )
    obsm = {k: v for k, v in obsm.items() if type(v) in obsm_types}
    varm = dict(
        array=np.random.random((N, 50)),
        sparse=sparse.random(N, 100, format="csr"),
        df=gen_typed_df(N, var_names),
    )
    varm = {k: v for k, v in varm.items() if type(v) in varm_types}
    layers = dict(
        array=np.random.random((M, N)), sparse=sparse.random(M, N, format="csr")
    )
    layers = {k: v for k, v in layers.items() if type(v) in layers_types}
    obsp = dict(
        array=np.random.random((M, M)), sparse=sparse.random(M, M, format="csr")
    )
    varp = dict(
        array=np.random.random((N, N)), sparse=sparse.random(N, N, format="csr")
    )
    uns = dict(
        O_recarray=gen_vstr_recarray(N, 5),
        nested=dict(
            scalar_str="str",
            scalar_int=42,
            scalar_float=3.0,
            nested_further=dict(array=np.arange(5)),
        ),
        # U_recarray=gen_vstr_recarray(N, 5, "U4")
    )
    adata = AnnData(
        X=X,
        obs=obs,
        var=var,
        obsm=obsm,
        varm=varm,
        layers=layers,
        obsp=obsp,
        varp=varp,
        dtype=X_dtype,
        uns=uns,
    )
    return adata


def array_bool_subset(index, min_size=2):
    b = np.zeros(len(index), dtype=bool)
    selected = np.random.choice(
        range(len(index)),
        size=np.random.randint(min_size, len(index), ()),
        replace=False,
    )
    b[selected] = True
    return b


def matrix_bool_subset(index, min_size=2):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PendingDeprecationWarning)
        indexer = np.matrix(
            array_bool_subset(index, min_size=min_size).reshape(len(index), 1)
        )
    return indexer


def spmatrix_bool_subset(index, min_size=2):
    return sparse.csr_matrix(
        array_bool_subset(index, min_size=min_size).reshape(len(index), 1)
    )


def array_subset(index, min_size=2):
    if len(index) < min_size:
        raise ValueError(
            f"min_size (={min_size}) must be smaller than len(index) (={len(index)}"
        )
    return np.random.choice(
        index, size=np.random.randint(min_size, len(index), ()), replace=False
    )


def array_int_subset(index, min_size=2):
    if len(index) < min_size:
        raise ValueError(
            f"min_size (={min_size}) must be smaller than len(index) (={len(index)}"
        )
    return np.random.choice(
        np.arange(len(index)),
        size=np.random.randint(min_size, len(index), ()),
        replace=False,
    )


def slice_subset(index, min_size=2):
    while True:
        points = np.random.choice(np.arange(len(index) + 1), size=2, replace=False)
        s = slice(*sorted(points))
        if len(range(*s.indices(len(index)))) >= min_size:
            break
    return s


def single_subset(index):
    return index[np.random.randint(0, len(index), size=())]


@pytest.fixture(
    params=[
        array_subset,
        slice_subset,
        single_subset,
        array_int_subset,
        array_bool_subset,
        matrix_bool_subset,
        spmatrix_bool_subset,
    ]
)
def subset_func(request):
    return request.param
