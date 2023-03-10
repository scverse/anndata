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
from itertools import product
import tempfile
from pathlib import Path
import pickle
import sys

from memory_profiler import memory_usage
import numpy as np
import pandas as pd
from scipy import sparse

from .utils import get_anndata_memsize, sedate, get_peak_mem, get_actualsize
#from . import datasets

import anndata


PBMC_3K_PATH = Path(__file__).parent / "data/pbmc3k_raw.h5ad"

PBMC_REDUCED_PATH = Path(__file__).parent / "10x_pbmc68k_reduced.h5ad"
BM_43K_CSR_PATH = Path(__file__).parent.parent / "datasets/BM2_43k-cells.h5ad"
BM_43K_CSC_PATH = Path(__file__).parent.parent / "datasets/BM2_43k-cells_CSC.h5ad"


# class ZarrReadSuite:
#     params = []
#     param_names = ["input_path"]

#     def setup(self, input_path):
#         self.filepath = input_path

#     def time_read_full(self, input_path):
#         anndata.read_zarr(self.filepath)

#     def peakmem_read_full(self, input_path):
#         anndata.read_zarr(self.filepath)

#     def mem_readfull_object(self, input_path):
#         return anndata.read_zarr(self.filepath)

#     def track_read_full_memratio(self, input_path):
#         mem_recording = memory_usage(
#             (sedate(anndata.read_zarr, 0.005), (self.filepath,)), interval=0.001
#         )
#         adata = anndata.read_zarr(self.filepath)
#         base_size = mem_recording[-1] - mem_recording[0]
#         print(np.max(mem_recording) - np.min(mem_recording))
#         print(base_size)
#         return (np.max(mem_recording) - np.min(mem_recording)) / base_size

#     def peakmem_read_backed(self, input_path):
#         anndata.read_zarr(self.filepath, backed="r")

#     def mem_read_backed_object(self, input_path):
#         return anndata.read_zarr(self.filepath, backed="r")


class H5ADInMemorySizeSuite:
    params = [PBMC_3K_PATH]
    param_names = ["input_path"]

    def setup(self, input_path):
        self.filepath = input_path

    def track_in_memory_size(self, input_path):
        adata = anndata.read_h5ad(self.filepath)
        adata_size = sys.getsizeof(adata)

        return adata_size

    def track_actual_in_memory_size(self, input_path):

        adata = anndata.read_h5ad(self.filepath)
        adata_size = get_actualsize(adata)

        return adata_size


class H5ADReadSuite:
    # params = [PBMC_REDUCED_PATH, PBMC_3K_PATH, BM_43K_CSR_PATH]
    params = [PBMC_3K_PATH]
    param_names = ["input_path"]

    def setup(self, input_path):
        self.filepath = input_path

    def time_read_full(self, input_path):
        anndata.read_h5ad(self.filepath)

    def peakmem_read_full(self, input_path):
        anndata.read_h5ad(self.filepath)

    def mem_readfull_object(self, input_path):
        return anndata.read_h5ad(self.filepath)

    def track_read_full_memratio(self, input_path):
        mem_recording = memory_usage(
            (sedate(anndata.read_h5ad, 0.005), (self.filepath,)), interval=0.001
        )
        # adata = anndata.read_h5ad(self.filepath)
        base_size = mem_recording[-1] - mem_recording[0]
        print(np.max(mem_recording) - np.min(mem_recording))
        print(base_size)
        return (np.max(mem_recording) - np.min(mem_recording)) / base_size

    def peakmem_read_backed(self, input_path):
        anndata.read_h5ad(self.filepath, backed="r")

    def mem_read_backed_object(self, input_path):
        return anndata.read_h5ad(self.filepath, backed="r")


class H5ADWriteSuite:
    # params = [PBMC_REDUCED_PATH, PBMC_3K_PATH, BM_43K_CSR_PATH]
    params = [PBMC_3K_PATH]
    param_names = ["input_path"]

    def setup(self, input_path):
        mem_recording, adata = memory_usage(
            (sedate(anndata.read_h5ad, 0.005), (input_path,)),
            retval=True,
            interval=0.001,
        )
        self.adata = adata
        self.base_size = mem_recording[-1] - mem_recording[0]
        self.tmpdir = tempfile.TemporaryDirectory()
        self.writepth = Path(self.tmpdir.name) / "out.h5ad"

    def teardown(self, input_path):
        self.tmpdir.cleanup()

    def time_write_full(self, input_path):
        self.adata.write_h5ad(self.writepth, compression=None)

    def peakmem_write_full(self, input_path):
        self.adata.write_h5ad(self.writepth)

    def track_peakmem_write_full(self, input_path):
        return get_peak_mem((sedate(self.adata.write_h5ad), (self.writepth,)))

    def time_write_compressed(self, input_path):
        self.adata.write_h5ad(self.writepth, compression="gzip")

    def peakmem_write_compressed(self, input_path):
        self.adata.write_h5ad(self.writepth, compression="gzip")

    def track_peakmem_write_compressed(self, input_path):
        return get_peak_mem(
            (sedate(self.adata.write_h5ad), (self.writepth,), {"compression": "gzip"})
        )


class H5ADBackedWriteSuite(H5ADWriteSuite):
    # params = [PBMC_REDUCED_PATH, PBMC_3K_PATH]
    params = [PBMC_3K_PATH]
    param_names = ["input_path"]

    def setup(self, input_path):
        mem_recording, adata = memory_usage(
            (sedate(anndata.read_h5ad, 0.005), (input_path,), {"backed": "r"}),
            retval=True,
            interval=0.001,
        )
        self.adata = adata
        self.base_size = mem_recording[-1] - mem_recording[0]
        self.tmpdir = tempfile.TemporaryDirectory()
        self.writepth = Path(self.tmpdir.name) / "out.h5ad"
