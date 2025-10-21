from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

import anndata as ad


class BackedHDF5:
    param_names = ("n_obs", "n_vars", "X_type")
    params = (
        (1000, 10000),
        (2000, 5000),
        ("dense", "sparse_csr"),
    )

    def setup(self, n_obs: int, n_vars: int, X_type: str):
        self.n_obs = n_obs
        self.n_vars = n_vars

        # Create count matrix
        if X_type == "dense":
            X = np.random.poisson(2, (n_obs, n_vars)).astype(np.float32)
        else:  # sparse_csr
            X = sparse.random(n_obs, n_vars, density=0.1, format="csr")
            X.data = np.random.poisson(5, X.data.shape).astype(np.float32)

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
                "gene_name": [f"gene_{i}" for i in range(n_vars)],
                "highly_variable": np.random.choice([True, False], n_vars),
            },
            index=[f"ENSG_{i:08d}" for i in range(n_vars)],
        )

        # Create AnnData object and save to HDF5
        adata = ad.AnnData(X=X, obs=obs, var=var)

        # Add some layers for more comprehensive testing
        adata.layers["raw"] = X.copy() if hasattr(X, "copy") else X
        adata.obsm["X_pca"] = np.random.randn(n_obs, 50)
        adata.varm["loadings"] = np.random.randn(n_vars, 50)

        # Create temporary file
        self.temp_path = Path(tempfile.mkdtemp()) / "test_data.h5ad"
        adata.write_h5ad(self.temp_path)

        # Open as backed
        self.adata_backed = ad.read_h5ad(self.temp_path, backed="r")

        # Prepare indices for duplicate index testing
        self.obs_idx_with_dupes = np.array([0, 1, 0, 2, 1] * (n_obs // 100 + 1))[
            : (n_obs // 10)
        ]
        self.var_idx_with_dupes = np.array([0, 1, 2, 0, 3] * (n_vars // 100 + 1))[
            : (n_vars // 10)
        ]
        self.obs_idx_no_dupes = np.arange(0, n_obs, 10)
        self.var_idx_no_dupes = np.arange(0, n_vars, 10)

    def time_slice_obs(self, *_):
        """Time slicing observations from backed HDF5"""
        self.adata_backed[0 : (self.n_obs // 2), :].X

    def peakmem_slice_obs(self, *_):
        """Peak memory for slicing observations from backed HDF5"""
        self.adata_backed[0 : (self.n_obs // 2), :].X

    def time_slice_vars(self, *_):
        """Time slicing variables from backed HDF5"""
        self.adata_backed[:, 0 : (self.n_vars // 2)].X

    def peakmem_slice_vars(self, *_):
        """Peak memory for slicing variables from backed HDF5"""
        self.adata_backed[:, 0 : (self.n_vars // 2)].X

    def time_fancy_index_no_dupes(self, *_):
        """Time fancy indexing without duplicates"""
        self.adata_backed[self.obs_idx_no_dupes, self.var_idx_no_dupes].X

    def peakmem_fancy_index_no_dupes(self, *_):
        """Peak memory for fancy indexing without duplicates"""
        self.adata_backed[self.obs_idx_no_dupes, self.var_idx_no_dupes].X

    def time_fancy_index_with_dupes_obs(self, *_):
        """Time fancy indexing with duplicate observation indices"""
        self.adata_backed[self.obs_idx_with_dupes, :].X

    def peakmem_fancy_index_with_dupes_obs(self, *_):
        """Peak memory for fancy indexing with duplicate observation indices"""
        self.adata_backed[self.obs_idx_with_dupes, :].X

    def time_fancy_index_with_dupes_vars(self, *_):
        """Time fancy indexing with duplicate variable indices"""
        self.adata_backed[:, self.var_idx_with_dupes].X

    def peakmem_fancy_index_with_dupes_vars(self, *_):
        """Peak memory for fancy indexing with duplicate variable indices"""
        self.adata_backed[:, self.var_idx_with_dupes].X

    def time_fancy_index_with_dupes_both(self, *_):
        """Time fancy indexing with duplicate indices in both dimensions"""
        self.adata_backed[self.obs_idx_with_dupes, self.var_idx_with_dupes].X

    def peakmem_fancy_index_with_dupes_both(self, *_):
        """Peak memory for fancy indexing with duplicate indices in both dimensions"""
        self.adata_backed[self.obs_idx_with_dupes, self.var_idx_with_dupes].X

    def time_boolean_indexing(self, *_):
        """Time boolean indexing"""
        mask = np.random.choice([True, False], self.n_obs, p=[0.1, 0.9])
        self.adata_backed[mask, :].X

    def peakmem_boolean_indexing(self, *_):
        """Peak memory for boolean indexing"""
        mask = np.random.choice([True, False], self.n_obs, p=[0.1, 0.9])
        self.adata_backed[mask, :].X

    def time_to_memory_subset(self, *_):
        """Time converting subset to memory"""
        subset = self.adata_backed[0 : (self.n_obs // 4), 0 : (self.n_vars // 4)]
        subset.to_memory()

    def peakmem_to_memory_subset(self, *_):
        """Peak memory for converting subset to memory"""
        subset = self.adata_backed[0 : (self.n_obs // 4), 0 : (self.n_vars // 4)]
        subset.to_memory()

    def time_access_layers(self, *_):
        """Time accessing layers from backed data"""
        self.adata_backed[0 : (self.n_obs // 10), :].layers["raw"]

    def peakmem_access_layers(self, *_):
        """Peak memory for accessing layers from backed data"""
        self.adata_backed[0 : (self.n_obs // 10), :].layers["raw"]

    def time_access_obsm(self, *_):
        """Time accessing obsm from backed data"""
        self.adata_backed[0 : (self.n_obs // 10), :].obsm["X_pca"]

    def peakmem_access_obsm(self, *_):
        """Peak memory for accessing obsm from backed data"""
        self.adata_backed[0 : (self.n_obs // 10), :].obsm["X_pca"]

    def teardown(self, *_):
        """Clean up temporary files"""
        if hasattr(self, "adata_backed"):
            self.adata_backed.file.close()
        if hasattr(self, "temp_path") and self.temp_path.exists():
            self.temp_path.unlink()
            self.temp_path.parent.rmdir()
