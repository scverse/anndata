import anndata as ad
import h5py
import pandas as pd
from scipy import sparse

from .utils import Dataset, download, save_sparse, load_sparse
from . import DATASETS_PATH

BM_PATH = DATASETS_PATH / "data/raw/ica_bone_marrow_h5.h5"


def fetch_bonemarrow():
    download(
        "https://s3.amazonaws.com/preview-ica-expression-data/ica_bone_marrow_h5.h5",
        BM_PATH,
    )


class ICA_BoneMarrow_full(Dataset):
    name = "ICA-BoneMarrow_378k-cells_34k-genes"
    dirpth = DATASETS_PATH / "data/ica_bone-marrow_full"

    @classmethod
    def load(cls):
        mtx = load_sparse(cls.dirpth / "mtx.h5")
        obs = pd.read_pickle(cls.dirpth / "obs.pkl.gz")
        var = pd.read_pickle(cls.dirpth / "var.pkl.gz")
        adata = ad.AnnData(X=mtx, obs=obs, var=var)
        return adata

    def setup(cls):
        if not BM_PATH.is_file():
            fetch_bonemarrow()
        if not cls.dirpth.is_dir():
            cls.dirpth.mkdir(parents=True, exist_ok=True)
        with h5py.File(BM_PATH, "r") as f:
            g = f["GRCh38"]
            save_sparse(
                cls.dirpth / "mtx.h5",
                sparse.csr_matrix(
                    (g["data"][:], g["indices"][:], g["indptr"][:]),
                    shape=g["shape"][:][::-1],
                ),
            )
            pd.DataFrame(
                {"gene_symbols": g["gene_names"][:].astype(str)},
                index=pd.Index(g["genes"][:]).astype(str),
            ).to_pickle(cls.dirpth / "var.pkl.gz")
            pd.DataFrame([], index=pd.Index(g["barcodes"][:].astype(str))).to_pickle(
                cls.dirpth / "obs.pkl.gz"
            )


class ICA_BoneMarrow_Donor1(Dataset):
    name = "ICA-BoneMarrow_48k-cells_34k-genes"
    dirpth = DATASETS_PATH / "data/ica_bone-marrow_donor1"

    @classmethod
    def load(cls):
        mtx = load_sparse(cls.dirpth / "mtx.h5")
        obs = pd.read_pickle(cls.dirpth / "obs.pkl.gz")
        var = pd.read_pickle(cls.dirpth / "var.pkl.gz")
        adata = ad.AnnData(X=mtx, obs=obs, var=var)
        return adata

    @classmethod
    def setup(cls):
        if not BM_PATH.is_file():
            fetch_bonemarrow()
        if not cls.dirpth.is_dir():
            cls.dirpth.mkdir(parents=True, exist_ok=True)
        with h5py.File(BM_PATH, "r") as f:
            g = f["GRCh38"]
            adata = ad.AnnData(
                X=sparse.csr_matrix(
                    (g["data"][:], g["indices"][:], g["indptr"][:]),
                    shape=g["shape"][:][::-1],
                ),
                obs=pd.DataFrame([], index=pd.Index(g["barcodes"][:]).astype(str)),
                var=pd.DataFrame(
                    {"gene_symbols": g["gene_names"][:].astype(str)},
                    index=pd.Index(g["genes"][:]).astype(str),
                ),
            )
            subset = adata[adata.obs_names.str.contains("BM1"), :].copy()
            save_sparse(cls.dirpth / "mtx.h5", subset.X)
            subset.var.to_pickle(cls.dirpth / "var.pkl.gz")
            subset.obs.to_pickle(cls.dirpth / "obs.pkl.gz")
