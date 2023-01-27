from pathlib import Path

import h5py
from scipy import sparse

###############################################################################
# Base classes for Datasets
###############################################################################


class MetaDataset(type):
    def __repr__(cls):
        return cls.name


class Dataset(object, metaclass=MetaDataset):
    @classmethod
    def is_available(cls):
        return cls.dirpth.is_dir()


###############################################################################
# IO
###############################################################################


def download(url: str, path: Path):
    from tqdm.auto import tqdm
    from urllib.request import urlretrieve

    path.parent.mkdir(parents=True, exist_ok=True)
    with tqdm(unit="B", unit_scale=True, miniters=1, desc=path.name) as t:

        def update_to(b=1, bsize=1, tsize=None):
            if tsize is not None:
                t.total = tsize
            t.update(b * bsize - t.n)

        urlretrieve(url, str(path), reporthook=update_to)


def save_sparse(pth, obj):
    with h5py.File(pth, "w") as f:
        f.create_dataset("data", data=obj.data, compression="gzip")
        f.create_dataset("indices", data=obj.indices, compression="gzip")
        f.create_dataset("indptr", data=obj.indptr, compression="gzip")
        f.attrs["shape"] = obj.shape
        f.attrs["format"] = obj.format


def load_sparse(pth):
    with h5py.File(pth, "r") as f:
        if f.attrs["format"] == "csr":
            mtx_class = sparse.csr_matrix
        elif f.attrs["format"] == "csc":
            mtx_class = sparse.csc_matrix
        else:
            raise NotImplementedError()
        mtx = mtx_class(
            (f["data"][:], f["indices"][:], f["indptr"][:]), shape=f.attrs["shape"]
        )
    return mtx
