from __future__ import annotations

import re
import warnings
from importlib.util import find_spec
from pathlib import Path

import h5py
import pytest
from packaging.version import Version

import anndata as ad
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, gen_adata


@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
def test_old_format_warning_thrown():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=r"Importing read_.* from `anndata` is deprecated"
        )
        import scanpy as sc

    def msg_re(entry: str) -> str:
        return re.escape(
            f"Moving element from .uns['neighbors'][{entry!r}] to .obsp[{entry!r}]."
        )

    pth = Path(sc.datasets.__file__).parent / "10x_pbmc68k_reduced.h5ad"
    with (
        pytest.warns(FutureWarning, match=msg_re("distances")),
        pytest.warns(FutureWarning, match=msg_re("connectivities")),
        pytest.warns(ad.OldFormatWarning),
    ):
        warnings.simplefilter("default", FutureWarning)
        ad.read_h5ad(pth)


def test_old_format_warning_not_thrown(tmp_path):
    pth = tmp_path / "current.h5ad"
    adata = gen_adata((20, 10), **GEN_ADATA_NO_XARRAY_ARGS)
    adata.write_h5ad(pth)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always", ad.OldFormatWarning)
        if Version(h5py.__version__) < Version("3.2"):
            # https://github.com/h5py/h5py/issues/1808
            warnings.filterwarnings(
                "ignore",
                r"Passing None into shape arguments as an alias for \(\) is deprecated\.",
                category=DeprecationWarning,
            )

        ad.read_h5ad(pth)

    if len(record) != 0:
        msg_content = "\n".join(
            [f"\t{w.category.__name__}('{w.message}')" for w in record]
        )
        pytest.fail(
            f"Warnings were thrown when they shouldn't be. Got:\n\n{msg_content}"
        )
