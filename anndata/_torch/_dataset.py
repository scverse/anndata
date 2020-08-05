import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.preprocessing import LabelEncoder

try:
    import torch
    from torch.utils.data import Dataset as DatasetTorch
    from torch.utils.data import DataLoader
except ImportError:
    print('could not load torch')
    torch = None
    DatasetTorch = object
    DataLoader = object


def data_loaders(data, label_key=None, batch_size=100, shuffle=False):
    """Generate `DataLoader` for each split from `AnnData`.

    Params
    ------
    data : `AnnData`, `Mapping[AnnData]`
        Needs to contain an `.obs` field named `split` defining the split.
    label_key : str
        Key for `.obs` defining the label.
    batch_size : int
        Batch size.
    shuffle : bool
        `DataLoader` `shuffle` parameter.

    Returns
    -------
    A `dict[DataLoader]` indexed by the names of splits.
    """
    datasets = {}  # torch data
    loaders = {}  # loaders

    if isinstance(data, dict):
        splits = list(data.keys())
    else:
        adata = data
        if 'split' not in adata.obs.columns:
            raise ValueError('Either pass dict with splits or an AnnData with obs column "split".')
        splits = list(adata.obs.split.unique())

    # check that we have training and test split
    if not ('train' in splits and 'test' in splits):
        raise ValueError('Need to have "train" and "test" in split.')
    # ensure train comes first for encoder
    splits.remove("train")
    splits.insert(0, "train")

    label_encoder = None  # is None for train
    for split in splits:
        if isinstance(data, dict):
            adata_split = data[split]
        else:
            if 'split' in adata.obs.columns:
                adata_split = adata[adata.obs.split == split].copy()  # make copy to avoid slow lazy indexing
        datasets[split] = Dataset(adata_split, label_key=label_key, label_encoder=label_encoder)
        # need to set the fitted label encoder so that it's used in validation,
        # test, holdout or whatever might come after in the loop
        if split == 'train':  # set once and never again
            label_encoder = datasets[split].label_encoder
        loaders[split] = DataLoader(
            dataset=datasets[split], batch_size=batch_size, shuffle=shuffle)

    return loaders


class Dataset(DatasetTorch):
    """torch `Dataset` from `AnnData`.

    Params
    ------
    adata : `AnnData`
        Needs to contain an `.obs` field named `split` defining the split.
    label_key : str
        Key for `.obs` defining the label.
    device : `torch.device`
    """
    def __init__(self, adata, label_key=None, label_encoder=None, device=None):
        # label key
        self._label_key = label_key
        # label encoder and label codes
        if label_key is not None:
            if label_encoder is None:
                label_encoder = LabelEncoder()
                codes = label_encoder.fit_transform(adata.obs[label_key].values)
            else:
                codes = label_encoder.transform(adata.obs[label_key].values)
            self._label_encoder = label_encoder
            # need to transform to int64 and 2d for torch
            self._label_codes = codes.astype('int64')[:, None]
        # adata
        if sparse.issparse(adata.X):
            print('warning: densifying data')
            adata = adata.copy()  # need to copy otherwise we'll change the outside object
            adata.X = adata.X.A  # this here shouldn't actually take memory?
        self._adata = adata
        # device
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    @property
    def adata(self):
        return self._adata

    @property
    def var_names(self):
        return self.adata.var_names

    @property
    def label_key(self):
        return self._label_key

    @property
    def label_codes(self):
        return self._label_codes

    @property
    def label_encoder(self):
        return self._label_encoder

    def __getitem__(self, index):
        X_tensor = torch.from_numpy(self.adata.X[index]).to(self.device)
        if self.label_key is not None:
            # need to transform to 2d
            label_tensor = torch.from_numpy(self.label_codes[index]).to(self.device)
            return X_tensor, label_tensor
        else:
            return X_tensor

    def __len__(self):
        return self.adata.n_obs

DatasetFromAnnData = Dataset  # backward compat
