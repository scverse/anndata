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
from types import MappingProxyType

import numpy as np
import pooch
from memory_profiler import memory_usage

# from . import datasets
import anndata

from .utils import get_actualsize, get_peak_mem, sedate

PBMC_3K_URL = "https://falexwolf.de/data/pbmc3k_raw.h5ad"

# PBMC_3K_PATH = Path(__file__).parent / "data/pbmc3k_raw.h5ad"
# PBMC_REDUCED_PATH = Path(__file__).parent / "10x_pbmc68k_reduced.h5ad"
# BM_43K_CSR_PATH = Path(__file__).parent.parent / "datasets/BM2_43k-cells.h5ad"
# BM_43K_CSC_PATH = Path(__file__).parent.parent / "datasets/BM2_43k-cells_CSC.h5ad"


# class ZarrReadSuite:
#     params = []
#     param_names = ["input_url"]

#     def setup(self, input_url):
#         self.filepath = pooch.retrieve(url=input_url, known_hash=None)

#     def time_read_full(self, input_url):
#         anndata.read_zarr(self.filepath)

#     def peakmem_read_full(self, input_url):
#         anndata.read_zarr(self.filepath)

#     def mem_readfull_object(self, input_url):
#         return anndata.read_zarr(self.filepath)

#     def track_read_full_memratio(self, input_url):
#         mem_recording = memory_usage(
#             (sedate(anndata.read_zarr, 0.005), (self.filepath,)), interval=0.001
#         )
#         adata = anndata.read_zarr(self.filepath)
#         base_size = mem_recording[-1] - mem_recording[0]
#         print(np.max(mem_recording) - np.min(mem_recording))
#         print(base_size)
#         return (np.max(mem_recording) - np.min(mem_recording)) / base_size

#     def peakmem_read_backed(self, input_url):
#         anndata.read_zarr(self.filepath, backed="r")

#     def mem_read_backed_object(self, input_url):
#         return anndata.read_zarr(self.filepath, backed="r")


class H5ADInMemorySizeSuite:
    _urls = MappingProxyType(dict(pbmc3k=PBMC_3K_URL))
    params = _urls.keys()
    param_names = ("input_data",)

    def setup(self, input_data: str):
        self.filepath = pooch.retrieve(url=self._urls[input_data], known_hash=None)

    def track_in_memory_size(self, *_):
        adata = anndata.read_h5ad(self.filepath)
        adata_size = sys.getsizeof(adata)

        return adata_size

    def track_actual_in_memory_size(self, *_):
        adata = anndata.read_h5ad(self.filepath)
        adata_size = get_actualsize(adata)

        return adata_size


class H5ADReadSuite:
    _urls = MappingProxyType(dict(pbmc3k=PBMC_3K_URL))
    params = _urls.keys()
    param_names = ("input_data",)

    def setup(self, input_data: str):
        self.filepath = pooch.retrieve(url=self._urls[input_data], known_hash=None)

    def time_read_full(self, *_):
        anndata.read_h5ad(self.filepath)

    def peakmem_read_full(self, *_):
        anndata.read_h5ad(self.filepath)

    def mem_readfull_object(self, *_):
        return anndata.read_h5ad(self.filepath)

    def track_read_full_memratio(self, *_):
        mem_recording = memory_usage(
            (sedate(anndata.read_h5ad, 0.005), (self.filepath,)), interval=0.001
        )
        # adata = anndata.read_h5ad(self.filepath)
        base_size = mem_recording[-1] - mem_recording[0]
        print(np.max(mem_recording) - np.min(mem_recording))
        print(base_size)
        return (np.max(mem_recording) - np.min(mem_recording)) / base_size

    def peakmem_read_backed(self, *_):
        anndata.read_h5ad(self.filepath, backed="r")

    # causes benchmarking to break from: https://github.com/pympler/pympler/issues/151
    # def mem_read_backed_object(self, *_):
    #     return anndata.read_h5ad(self.filepath, backed="r")


class H5ADWriteSuite:
    _urls = MappingProxyType(dict(pbmc3k=PBMC_3K_URL))
    params = _urls.keys()
    param_names = ("input_data",)

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
        self.adata.write_h5ad(self.writepth, compression=None)

    def peakmem_write_full(self, *_):
        self.adata.write_h5ad(self.writepth)

    def track_peakmem_write_full(self, *_):
        return get_peak_mem((sedate(self.adata.write_h5ad), (self.writepth,)))

    def time_write_compressed(self, *_):
        self.adata.write_h5ad(self.writepth, compression="gzip")

    def peakmem_write_compressed(self, *_):
        self.adata.write_h5ad(self.writepth, compression="gzip")

    def track_peakmem_write_compressed(self, *_):
        return get_peak_mem(
            (sedate(self.adata.write_h5ad), (self.writepth,), {"compression": "gzip"})
        )


class H5ADBackedWriteSuite(H5ADWriteSuite):
    _urls = MappingProxyType(dict(pbmc3k=PBMC_3K_URL))
    params = _urls.keys()
    param_names = ("input_data",)

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
