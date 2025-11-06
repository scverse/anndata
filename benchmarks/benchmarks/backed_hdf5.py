from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

import anndata as ad

file_paths = {"sparse": "adata_sparse.h5ad"}


class BackedHDF5Indexing:
    param_names = ("arr_type",)
    params = ("sparse",)

    def setup_cache(self) -> None:
        x_sparse = sparse.random(
            10000,
            50000,
            density=0.01,
            format="csr",
            random_state=np.random.default_rng(42),
        )
        for X, arr_type in [
            (x_sparse, "sparse"),
        ]:
            n_obs, n_var = X.shape

            # Create obs and var dataframes
            obs = pd.DataFrame(
                {
                    "cell_type": pd.Categorical(
                        np.random.choice(["TypeA", "TypeB", "TypeC"], n_obs)
                    ),
                    "total_counts": np.random.randint(1000, 5000, n_obs),
                },
                index=[f"cell_{i}" for i in range(n_obs)],
            )

            var = pd.DataFrame(
                {
                    "gene_name": [f"gene_{i}" for i in range(n_var)],
                },
                index=[f"ENSG_{i:08d}" for i in range(n_var)],
            )

            # Create AnnData object and save to HDF5
            adata = ad.AnnData(X=X, obs=obs, var=var)

            # Create temporary file
            adata.write_h5ad(file_paths[arr_type])

    def setup(self, arr_type):
        # Open as backed
        self.adata_backed = ad.read_h5ad(file_paths[arr_type], backed="r")
        self.n_obs, self.n_var = self.adata_backed.shape
        # Prepare indices for duplicate index testing
        self.obs_idx_with_dupes = np.array([0, 1, 0, 2, 1] * (self.n_obs // 100 + 1))[
            : (self.n_obs // 10)
        ]
        self.var_idx_with_dupes = np.array([0, 1, 2, 0, 3] * (self.n_var // 100 + 1))[
            : (self.n_var // 10)
        ]
        self.obs_idx_no_dupes = np.arange(0, self.n_obs, 10)
        self.var_idx_no_dupes = np.arange(0, self.n_var, 10)

    def time_slice_obs(self, *_):
        """Time slicing observations from backed HDF5"""
        self.adata_backed[0 : (self.n_obs // 2), :]

    def time_slice_obs_to_memory(self, *_):
        """Time slicing observations from backed HDF5"""
        self.adata_backed[0 : (self.n_obs // 2), :].to_memory()

    def peakmem_slice_obs(self, *_):
        """Peak memory for slicing observations from backed HDF5"""
        self.adata_backed[0 : (self.n_obs // 2), :]

    def time_fancy_index_no_dupes(self, *_):
        """Time fancy indexing without duplicates"""
        self.adata_backed[self.obs_idx_no_dupes, self.var_idx_no_dupes]

    def peakmem_fancy_index_no_dupes(self, *_):
        """Peak memory for fancy indexing without duplicates"""
        self.adata_backed[self.obs_idx_no_dupes, self.var_idx_no_dupes]

    def time_fancy_index_no_dupes_to_memory(self, *_):
        """Time fancy indexing without duplicates"""
        self.adata_backed[self.obs_idx_no_dupes, self.var_idx_no_dupes].to_memory()

    def time_index_with_dupes_obs(self, *_):
        """Time fancy indexing with duplicate observation indices"""
        self.adata_backed[self.obs_idx_with_dupes, :]

    def peakmem_index_with_dupes_obs(self, *_):
        """Peak memory for fancy indexing with duplicate observation indices"""
        self.adata_backed[self.obs_idx_with_dupes, :]

    def time_to_memory_subset(self, *_):
        """Time converting subset to memory"""
        subset = self.adata_backed[0 : (self.n_obs // 4), 0 : (self.n_var // 4)]
        subset.to_memory()

    def peakmem_to_memory_subset(self, *_):
        """Peak memory for converting subset to memory"""
        subset = self.adata_backed[0 : (self.n_obs // 4), 0 : (self.n_var // 4)]
        subset.to_memory()

    def teardown(self, *_):
        """Clean up temporary files"""
        if hasattr(self, "adata_backed"):
            self.adata_backed.file.close()
