class WriteWarning(UserWarning):
    pass


class OldFormatWarning(PendingDeprecationWarning):
    pass


from .read import (
    read_csv,
    read_excel,
    read_umi_tools,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_zarr,
    read_h5ad,
)
from .write import write_csvs, write_loom, _write_h5ad, write_zarr
from . import h5ad
