"""Annotated multivariate observation data."""

from ._metadata import __version__, within_flit

if not within_flit():
    del within_flit
    from ._core.anndata import AnnData
    from ._core.merge import concat
    from ._core.raw import Raw
    from ._io import (
        read_h5ad,
        read_loom,
        read_hdf,
        read_excel,
        read_umi_tools,
        read_csv,
        read_text,
        read_mtx,
        read_zarr,
        read_10x_h5,
        read_10x_mtx,
    )
    from ._warnings import (
        OldFormatWarning,
        WriteWarning,
        ImplicitModificationWarning,
        ExperimentalFeatureWarning,
    )

    # backwards compat / shortcut for default format
    from ._io import read_h5ad as read
