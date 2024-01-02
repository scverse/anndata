from __future__ import annotations

import warnings
from importlib.util import find_spec
from pathlib import Path

import pytest

import anndata as ad
from anndata.tests.helpers import gen_adata


@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
def test_old_format_warning_thrown():
    import scanpy as sc

    with pytest.warns(ad._warnings.OldFormatWarning):
        pth = Path(sc.datasets.__file__).parent / "10x_pbmc68k_reduced.h5ad"
        ad.read_h5ad(pth)


def test_old_format_warning_not_thrown(tmp_path):
    pth = tmp_path / "current.h5ad"
    adata = gen_adata((20, 10))
    adata.write_h5ad(pth)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always", ad.OldFormatWarning)

        ad.read_h5ad(pth)

    if len(record) != 0:
        msg_content = "\n".join(
            [f"\t{w.category.__name__}('{w.message}')" for w in record]
        )
        pytest.fail(
            f"Warnings were thrown when they shouldn't be. Got:\n\n{msg_content}"
        )
