import warnings

import pytest

import anndata as ad
from anndata.tests.helpers import gen_adata


def test_old_format_warning_thrown():
    import scanpy as sc

    with pytest.warns(ad._io.OldFormatWarning):
        sc.datasets.pbmc68k_reduced()


def test_old_format_warning_not_thrown(tmp_path):
    pth = tmp_path / "current.h5ad"
    adata = gen_adata((20, 10))
    adata.write_h5ad(pth)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always", ad._io.OldFormatWarning)

        ad.read_h5ad(pth)

    if len(record) != 0:
        msg_content = "\n".join(
            [f"\t{w.category.__name__}('{w.message}')" for w in record]
        )
        pytest.fail(
            f"Warnings were thrown when they shouldn't be. Got:\n\n{msg_content}"
        )
