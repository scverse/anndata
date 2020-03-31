import numpy as np
from numpy import ma
import pytest
from scipy import sparse


from anndata import AnnData, Raw
from anndata.tests.helpers import assert_equal


def test_concatenate_dense():
    # dense data
    X1 = np.array([[1, 2, 3], [4, 5, 6]])
    X2 = np.array([[1, 2, 3], [4, 5, 6]])
    X3 = np.array([[1, 2, 3], [4, 5, 6]])

    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c"], annoA=[0, 1, 2]),
        obsm=dict(X_1=X1, X_2=X2, X_3=X3),
        layers=dict(Xs=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        obsm=dict(X_1=X1, X_2=X2, X_3=X3),
        layers={"Xs": X2},
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s1", "s2"], anno2=["d3", "d4"]),
        dict(var_names=["d", "c", "b"], annoB=[0, 1, 2]),
        obsm=dict(X_1=X1, X_2=X2),
        layers=dict(Xs=X3),
    )

    # inner join
    adata = adata1.concatenate(adata2, adata3)
    X_combined = [[2, 3], [5, 6], [3, 2], [6, 5], [3, 2], [6, 5]]
    assert adata.X.astype(int).tolist() == X_combined
    assert adata.layers["Xs"].astype(int).tolist() == X_combined
    assert adata.obs_keys() == ["anno1", "anno2", "batch"]
    assert adata.var_keys() == ["annoA-0", "annoA-1", "annoB-2"]
    assert adata.var.values.tolist() == [[1, 2, 2], [2, 1, 1]]
    assert adata.obsm_keys() == ["X_1", "X_2"]
    assert adata.obsm["X_1"].tolist() == np.concatenate([X1, X1, X1]).tolist()

    # with batch_key and batch_categories
    adata = adata1.concatenate(adata2, adata3, batch_key="batch1")
    assert adata.obs_keys() == ["anno1", "anno2", "batch1"]
    adata = adata1.concatenate(adata2, adata3, batch_categories=["a1", "a2", "a3"])
    assert adata.obs["batch"].cat.categories.tolist() == ["a1", "a2", "a3"]
    assert adata.var_names.tolist() == ["b", "c"]

    # outer join
    adata = adata1.concatenate(adata2, adata3, join="outer")

    Xma = ma.masked_invalid(adata.X)
    Xma_ref = ma.masked_invalid(
        np.array(
            [
                [1.0, 2.0, 3.0, np.nan],
                [4.0, 5.0, 6.0, np.nan],
                [np.nan, 3.0, 2.0, 1.0],
                [np.nan, 6.0, 5.0, 4.0],
                [np.nan, 3.0, 2.0, 1.0],
                [np.nan, 6.0, 5.0, 4.0],
            ]
        )
    )
    assert np.array_equal(Xma.mask, Xma_ref.mask)
    assert np.allclose(Xma.compressed(), Xma_ref.compressed())
    var_ma = ma.masked_invalid(adata.var.values.tolist())
    var_ma_ref = ma.masked_invalid(
        np.array(
            [
                [0.0, np.nan, np.nan],
                [1.0, 2.0, 2.0],
                [2.0, 1.0, 1.0],
                [np.nan, 0.0, 0.0],
            ]
        )
    )
    assert np.array_equal(var_ma.mask, var_ma_ref.mask)
    assert np.allclose(var_ma.compressed(), var_ma_ref.compressed())


def test_concatenate_dense_duplicates():
    X1 = np.array([[1, 2, 3], [4, 5, 6]])
    X2 = np.array([[1, 2, 3], [4, 5, 6]])
    X3 = np.array([[1, 2, 3], [4, 5, 6]])

    # inner join duplicates
    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(
            var_names=["a", "b", "c"],
            annoA=[0, 1, 2],
            annoB=[1.1, 1.0, 2.0],
            annoC=[1.1, 1.0, 2.0],
            annoD=[2.1, 2.0, 3.0],
        ),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(
            var_names=["a", "b", "c"],
            annoA=[0, 1, 2],
            annoB=[1.1, 1.0, 2.0],
            annoC=[1.1, 1.0, 2.0],
            annoD=[2.1, 2.0, 3.0],
        ),
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s1", "s2"], anno2=["d3", "d4"]),
        dict(
            var_names=["a", "b", "c"],
            annoA=[0, 1, 2],
            annoB=[1.1, 1.0, 2.0],
            annoD=[2.1, 2.0, 3.1],
        ),
    )

    adata = adata1.concatenate(adata2, adata3)
    assert adata.var_keys() == [
        "annoA",
        "annoB",
        "annoC-0",
        "annoD-0",
        "annoC-1",
        "annoD-1",
        "annoD-2",
    ]


def test_concatenate_sparse():
    # sparse data
    from scipy.sparse import csr_matrix

    X1 = csr_matrix([[0, 2, 3], [0, 5, 6]])
    X2 = csr_matrix([[0, 2, 3], [0, 5, 6]])
    X3 = csr_matrix([[1, 2, 0], [0, 5, 6]])

    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c"]),
        layers=dict(Xs=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(var_names=["d", "c", "b"]),
        layers=dict(Xs=X2),
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s5", "s6"], anno2=["d3", "d4"]),
        dict(var_names=["d", "c", "b"]),
        layers=dict(Xs=X3),
    )

    # inner join
    adata = adata1.concatenate(adata2, adata3)
    X_combined = [[2, 3], [5, 6], [3, 2], [6, 5], [0, 2], [6, 5]]
    assert adata.X.toarray().astype(int).tolist() == X_combined
    assert adata.layers["Xs"].toarray().astype(int).tolist() == X_combined

    # outer join
    adata = adata1.concatenate(adata2, adata3, join="outer")
    assert adata.X.toarray().tolist() == [
        [0.0, 2.0, 3.0, 0.0],
        [0.0, 5.0, 6.0, 0.0],
        [0.0, 3.0, 2.0, 0.0],
        [0.0, 6.0, 5.0, 0.0],
        [0.0, 0.0, 2.0, 1.0],
        [0.0, 6.0, 5.0, 0.0],
    ]


def test_concatenate_mixed():
    X1 = sparse.csr_matrix(np.array([[1, 2, 0], [4, 0, 6], [0, 0, 9]]))
    X2 = sparse.csr_matrix(np.array([[0, 2, 3], [4, 0, 0], [7, 0, 9]]))
    X3 = sparse.csr_matrix(np.array([[1, 0, 3], [0, 0, 6], [0, 8, 0]]))
    X4 = np.array([[0, 2, 3], [4, 0, 0], [7, 0, 9]])
    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2", "s3"], anno1=["c1", "c2", "c3"]),
        dict(var_names=["a", "b", "c"], annoA=[0, 1, 2]),
        layers=dict(counts=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s4", "s5", "s6"], anno1=["c3", "c4", "c5"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        layers=dict(counts=X4),  # sic
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s7", "s8", "s9"], anno2=["d3", "d4", "d5"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 2, 3], annoB=[0, 1, 2]),
        layers=dict(counts=X3),
    )
    adata4 = AnnData(
        X4,
        dict(obs_names=["s4", "s5", "s6"], anno1=["c3", "c4", "c5"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        layers=dict(counts=X2),  # sic
    )

    adata_all = AnnData.concatenate(adata1, adata2, adata3, adata4)
    assert isinstance(adata_all.X, sparse.csr_matrix)
    assert isinstance(adata_all.layers["counts"], sparse.csr_matrix)


def test_concatenate_with_raw():
    # dense data
    X1 = np.array([[1, 2, 3], [4, 5, 6]])
    X2 = np.array([[1, 2, 3], [4, 5, 6]])
    X3 = np.array([[1, 2, 3], [4, 5, 6]])

    X4 = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])

    adata1 = AnnData(
        X1,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c"], annoA=[0, 1, 2]),
        layers=dict(Xs=X1),
    )
    adata2 = AnnData(
        X2,
        dict(obs_names=["s3", "s4"], anno1=["c3", "c4"]),
        dict(var_names=["d", "c", "b"], annoA=[0, 1, 2]),
        layers=dict(Xs=X2),
    )
    adata3 = AnnData(
        X3,
        dict(obs_names=["s1", "s2"], anno2=["d3", "d4"]),
        dict(var_names=["d", "c", "b"], annoB=[0, 1, 2]),
        layers=dict(Xs=X3),
    )

    adata4 = AnnData(
        X4,
        dict(obs_names=["s1", "s2"], anno1=["c1", "c2"]),
        dict(var_names=["a", "b", "c", "z"], annoA=[0, 1, 2, 3]),
        layers=dict(Xs=X4),
    )

    adata1.raw = adata1
    adata2.raw = adata2
    adata3.raw = adata3

    adata_all = AnnData.concatenate(adata1, adata2, adata3)
    assert isinstance(adata_all.raw, Raw)
    assert set(adata_all.raw.var_names) == {"b", "c"}
    assert_equal(adata_all.raw.to_adata().obs, adata_all.obs)
    assert np.array_equal(adata_all.raw.X, adata_all.X)

    adata_all = AnnData.concatenate(adata1, adata2, adata3, join="outer")
    assert isinstance(adata_all.raw, Raw)
    assert set(adata_all.raw.var_names) == set("abcd")
    assert_equal(adata_all.raw.to_adata().obs, adata_all.obs)
    assert np.array_equal(np.nan_to_num(adata_all.raw.X), np.nan_to_num(adata_all.X))

    adata3.raw = adata4
    adata_all = AnnData.concatenate(adata1, adata2, adata3, join="outer")
    assert isinstance(adata_all.raw, Raw)
    assert set(adata_all.raw.var_names) == set("abcdz")
    assert set(adata_all.var_names) == set("abcd")
    assert not np.array_equal(
        np.nan_to_num(adata_all.raw.X), np.nan_to_num(adata_all.X)
    )

    del adata3.raw
    with pytest.warns(
        UserWarning,
        match=(
            "Only some adata objects have `.raw` attribute, "
            "not concatenating `.raw` attributes."
        ),
    ):
        adata_all = AnnData.concatenate(adata1, adata2, adata3)
    assert adata_all.raw is None

    del adata1.raw
    del adata2.raw
    assert all(_adata.raw is None for _adata in (adata1, adata2, adata3))
    adata_all = AnnData.concatenate(adata1, adata2, adata3)
    assert adata_all.raw is None
