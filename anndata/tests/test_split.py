import pandas as pd
import scanpy as sc
from anndata import split


def test_split():
    adata = sc.datasets.pbmc68k_reduced()

    adatas = split(adata, "bulk_labels")
    for g, ad in adatas.items():
        assert ad.obs["bulk_labels"].cat.categories.tolist() == [g]
        assert ad.is_view

    groups = ["CD14+ Monocyte", "CD34+"]
    adatas = split(adata, "bulk_labels", groups, copy=True)
    assert list(adatas.keys()) == groups
    for g in groups:
        assert adatas[g].obs["bulk_labels"].cat.categories.tolist() == [g]
        assert not adatas[g].is_view

    adatas = split(
        adata,
        "bulk_labels",
        dict(some=groups),
        others_key="others",
    )
    assert list(adatas.keys()) == ["some", "others"]
    assert all(adatas["some"].obs["bulk_labels"].isin(groups))
    mask = ~adata.obs["bulk_labels"].isin(groups)
    assert adatas["others"].obs_names.equals(adata[mask].obs_names)

    var_cats = pd.cut(adata.var.n_counts, 4).cat.rename_categories(str)
    adatas = split(adata, var_cats, axis=1)
    assert set(adatas.keys()) == set(var_cats.cat.categories)
    assert sum(a.n_vars for a in adatas.values()) == adata.n_vars
