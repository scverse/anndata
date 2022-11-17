import anndata as ad
import h5py
import pandas as pd
from scipy import sparse

from .utils import Dataset, download, save_sparse, load_sparse
from . import DATASETS_PATH

PBMC_10x_NEXTGEM_PATH = (
    DATASETS_PATH / "data/raw/5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5"
)


def fetch_pbmc():
    download(
        "http://cf.10xgenomics.com/samples/cell-exp/3.1.0/5k_pbmc_protein_v3_nextgem/5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5",
        PBMC_10x_NEXTGEM_PATH,
    )


class PBMC_10X_NextGEM(Dataset):
    name = "10X-PBMC_5k-cells_34k-genes"
    dirpth = DATASETS_PATH / "data/10X_PBMC_NextGEM"

    @classmethod
    def load(cls):
        mtx = load_sparse(cls.dirpth / "mtx.h5")
        obs = pd.read_pickle(cls.dirpth / "obs.pkl.gz")
        var = pd.read_pickle(cls.dirpth / "var.pkl.gz")
        adata = ad.AnnData(X=mtx, obs=obs, var=var)
        return adata

    @classmethod
    def setup(cls):
        if not PBMC_10x_NEXTGEM_PATH.is_file():
            fetch_pbmc()
        if not cls.dirpth.is_dir():
            cls.dirpth.mkdir(parents=True, exist_ok=True)
        with h5py.File(PBMC_10x_NEXTGEM_PATH, "r") as f:
            g = f["matrix"]
            features = g["features"]
            adata = ad.AnnData(
                X=sparse.csr_matrix(
                    (g["data"][:], g["indices"][:], g["indptr"][:]),
                    shape=g["shape"][:][::-1],
                ),
                obs=pd.DataFrame([], index=pd.Index(g["barcodes"][:]).astype(str)),
                var=pd.DataFrame(
                    {
                        "gene_symbols": features["name"][:].astype(str),
                        "feature_type": pd.Categorical(
                            features["feature_type"][:].astype(str)
                        ),
                    },
                    index=pd.Index(features["id"][:].astype(str)),
                ),
            )

        subset = adata[:, adata.var["feature_type"] == "Gene Expression"].copy()
        subset.var.drop(columns="feature_type", inplace=True)
        save_sparse(cls.dirpth / "mtx.h5", subset.X)
        subset.var.to_pickle(cls.dirpth / "var.pkl.gz")
        subset.obs.to_pickle(cls.dirpth / "obs.pkl.gz")
