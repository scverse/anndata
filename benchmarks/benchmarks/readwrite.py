"""
This module will benchmark io of AnnData objects

Things to test:

* Read time, write time
* Peak memory during io
* File sizes

Parameterized by:

* What method is being used
* What data is being included
* Size of data being used

Also interesting:

* io for views
* io for backed objects
* Reading dense as sparse, writing sparse as dense
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pooch
from memory_profiler import memory_usage

# from . import datasets
import anndata

from .utils import get_actualsize, get_peak_mem, sedate

if TYPE_CHECKING:
    from collections.abc import Callable

PBMC_3K_URL = "https://falexwolf.de/data/pbmc3k_raw.h5ad"

# PBMC_3K_PATH = Path(__file__).parent / "data/pbmc3k_raw.h5ad"
# PBMC_REDUCED_PATH = Path(__file__).parent / "10x_pbmc68k_reduced.h5ad"
# BM_43K_CSR_PATH = Path(__file__).parent.parent / "datasets/BM2_43k-cells.h5ad"
# BM_43K_CSC_PATH = Path(__file__).parent.parent / "datasets/BM2_43k-cells_CSC.h5ad"


class TestSuite:
    _urls = dict(pbmc3k=PBMC_3K_URL)
    params = _urls.keys()
    param_names = ["input_data"]
    filepath: Path
    read_func: Callable[[Path | str], anndata.AnnData]

    def setup(self, input_data: str):
        self.filepath = Path(
            pooch.retrieve(url=self._urls[input_data], known_hash=None)
        )


class ZarrMixin(TestSuite):
    def setup(self, input_data: str):
        super().setup(input_data)
        zarr_path = self.filepath.with_suffix(".zarr")
        anndata.read_h5ad(self.filepath).write_zarr(zarr_path)
        self.filepath = zarr_path

    @property
    def read_func(self):
        return anndata.read_zarr


class H5ADInMemorySizeSuite(TestSuite):
    @property
    def read_func(self):
        return anndata.read_h5ad

    def track_in_memory_size(self, *_):
        adata = self.read_func(self.filepath)
        adata_size = sys.getsizeof(adata)

        return adata_size

    def track_actual_in_memory_size(self, *_):
        adata = self.read_func(self.filepath)
        adata_size = get_actualsize(adata)

        return adata_size


class ZarrInMemorySizeSuite(ZarrMixin, H5ADInMemorySizeSuite):
    @property
    def read_func(self):
        return anndata.read_zarr


class H5ADReadSuite(TestSuite):
    @property
    def read_func(self):
        return anndata.read_h5ad

    def time_read_full(self, *_):
        self.read_func(self.filepath)

    def peakmem_read_full(self, *_):
        self.read_func(self.filepath)

    def mem_readfull_object(self, *_):
        return self.read_func(self.filepath)

    def track_read_full_memratio(self, *_):
        mem_recording = memory_usage(
            (sedate(self.read_func, 0.005), (self.filepath,)), interval=0.001
        )
        # adata = self.read_func(self.filepath)
        base_size = mem_recording[-1] - mem_recording[0]
        print(np.max(mem_recording) - np.min(mem_recording))
        print(base_size)
        return (np.max(mem_recording) - np.min(mem_recording)) / base_size

    # causes benchmarking to break from: https://github.com/pympler/pympler/issues/151
    # def mem_read_backed_object(self, *_):
    #     return self.read_func(self.filepath, backed="r")


class BackedH5ADSuite(TestSuite):
    def peakmem_read_backed(self, *_):
        anndata.read_h5ad(self.filepath, backed="r")


class ZarrReadSuite(ZarrMixin, H5ADReadSuite):
    @property
    def read_func(self):
        return anndata.read_zarr


class H5ADWriteSuite:
    _urls = dict(pbmc3k=PBMC_3K_URL)
    params = _urls.keys()
    param_names = ["input_data"]
    write_func_str = "write_h5ad"

    def setup(self, input_data: str):
        mem_recording, adata = memory_usage(
            (
                sedate(anndata.read_h5ad, 0.005),
                (pooch.retrieve(self._urls[input_data], known_hash=None),),
            ),
            retval=True,
            interval=0.001,
        )
        self.adata = adata
        self.base_size = mem_recording[-1] - mem_recording[0]
        self.tmpdir = tempfile.TemporaryDirectory()
        self.writepth = Path(self.tmpdir.name) / "out.h5ad"

    def teardown(self, *_):
        self.tmpdir.cleanup()

    def time_write_full(self, *_):
        getattr(self.adata, self.write_func_str)(self.writepth, compression=None)

    def peakmem_write_full(self, *_):
        getattr(self.adata, self.write_func_str)(self.writepth)

    def track_peakmem_write_full(self, *_):
        return get_peak_mem(
            (sedate(getattr(self.adata, self.write_func_str)), (self.writepth,))
        )

    def time_write_compressed(self, *_):
        getattr(self.adata, self.write_func_str)(self.writepth, compression="gzip")

    def peakmem_write_compressed(self, *_):
        getattr(self.adata, self.write_func_str)(self.writepth, compression="gzip")

    def track_peakmem_write_compressed(self, *_):
        return get_peak_mem(
            (
                sedate(getattr(self.adata, self.write_func_str)),
                (self.writepth,),
                {"compression": "gzip"},
            )
        )


class ZarrWriteSizeSuite(H5ADWriteSuite):
    write_func_str = "write_zarr"

    def setup(self, input_data: str):
        h5_path = Path(pooch.retrieve(self._urls[input_data], known_hash=None))
        zarr_path = h5_path.with_suffix(".zarr")
        anndata.read_h5ad(h5_path).write_zarr(zarr_path)

        mem_recording, adata = memory_usage(
            (
                sedate(anndata.read_zarr, 0.005),
                (zarr_path,),
            ),
            retval=True,
            interval=0.001,
        )
        self.adata = adata
        self.base_size = mem_recording[-1] - mem_recording[0]
        self.tmpdir = tempfile.TemporaryDirectory()
        self.writepth = Path(self.tmpdir.name) / "out.zarr"

    def track_peakmem_write_compressed(self, *_):
        return get_peak_mem(
            (
                sedate(anndata.write_zarr, self.adata),
                (self.writepth,),
                {"compression": "gzip"},
            )
        )


class BackedH5ADWriteSuite(H5ADWriteSuite):
    def setup(self, input_data):
        mem_recording, adata = memory_usage(
            (
                sedate(anndata.read_h5ad, 0.005),
                (pooch.retrieve(self._urls[input_data], known_hash=None),),
                {"backed": "r"},
            ),
            retval=True,
            interval=0.001,
        )
        self.adata = adata
        self.base_size = mem_recording[-1] - mem_recording[0]
        self.tmpdir = tempfile.TemporaryDirectory()
        self.writepth = Path(self.tmpdir.name) / "out.h5ad"
