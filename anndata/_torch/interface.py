from sklearn.preprocessing import LabelEncoder
from scipy.sparse import issparse
from .._core.anndata import AnnData

try:
    import torch
    from torch.utils.data import Sampler, Dataset, DataLoader
except ImportError:
    raise ImportError("Please install pytorch")

# Custom sampler to get proper batches instead of joined separate indices
class BatchIndexSampler(Sampler):
    def __init__(self, n_obs, batch_size, shuffle=False, drop_last=False):
        self.n_obs = n_obs
        self.batch_size = batch_size if batch_size < n_obs else n_obs
        self.shuffle = shuffle
        self.drop_last = drop_last

    def __iter__(self):
        if self.shuffle:
            indices = torch.randperm(self.n_obs).tolist()
        else:
            indices = list(range(self.n_obs))

        for i in range(0, self.n_obs, self.batch_size):
            batch = indices[i : min(i + self.batch_size, self.n_obs)]

            # only happens if the last batch is smaller then batch_size
            if len(batch) < self.batch_size and self.drop_last:
                continue

            yield batch

    def __len__(self):
        from math import ceil

        if self.drop_last:
            length = self.n_obs // self.batch_size
        else:
            length = ceil(self.n_obs / self.batch_size)

        return length


class AnnDataSet(Dataset):
    def __init__(
        self,
        adata,
        obsm=None,
        layer=None,
        label_key=None,
        label_encoder=None,
        device=None,
        densify=False,
    ):
        self._adata = adata
        self._obsm = obsm
        self._layer = layer

        if obsm is None and layer is None:
            self._data = self._adata.X
        elif obsm is None and layer is not None:
            self._data = adata.layers[layer]
        elif obsm is not None and layer is None:
            self._data = adata.obsm[obsm]
        else:
            raise ValueError("You can't specify both obsm and layer.")

        self._sparse = issparse(self._data)

        if self._sparse and densify:
            self._data = self._data.toarray()
            self._sparse = False

        # copy from Alex's code
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
            self._label_codes = codes.astype("int64")[:, None]

        self.device = device

    def __getitem__(self, index):
        if self._sparse:
            data_tensor = torch.from_numpy(self._data[index].toarray())
        else:
            data_tensor = torch.from_numpy(self._data[index])

        label_tensor = None
        if self._label_key is not None:
            label_tensor = torch.from_numpy(self._label_codes[index])

        if self.device is not None:
            data_tensor = data_tensor.to(device)
            label_tensor = label_tensor.to(device) if label_tensor is not None else None

        if label_tensor is not None:
            return data_tensor, label_tensor
        else:
            return data_tensor

    def __len__(self):
        return self._adata.n_obs


# AnnDataLoader has the same arguments as DataLoader, but uses BatchIndexSampler by default
class AnnDataLoader(DataLoader):
    def __init__(self, adata, batch_size=1, shuffle=False, **kwargs):

        if isinstance(adata, AnnData):
            obsm = kwargs.pop("obsm", None)
            layer = kwargs.pop("layer", None)
            label_key = kwargs.pop("label_key", None)
            label_encoder = kwargs.pop("label_encoder", None)
            device = kwargs.pop("device", None)
            densify = kwargs.pop("densify", False)

            dataset = AnnDataSet(
                adata, obsm, layer, label_key, label_encoder, device, densify
            )

        elif isinstance(adata, AnnDataSet):
            dataset = adata
        else:
            raise ValueError("adata should be of type AnnData or AnnDataSet.")

        has_sampler = "sampler" in kwargs
        has_batch_sampler = "batch_sampler" in kwargs
        has_worker_init_fn = "worker_init_fn" in kwargs

        if (
            batch_size is not None
            and batch_size > 1
            and not has_sampler
            and not has_batch_sampler
            and not has_worker_init_fn
        ):
            drop_last = kwargs.pop("drop_last", False)
            default_sampler = BatchIndexSampler(
                len(dataset), batch_size, shuffle, drop_last
            )

            super().__init__(
                dataset, batch_size=None, sampler=default_sampler, **kwargs
            )
        else:
            super().__init__(dataset, batch_size=batch_size, shuffle=shuffle, **kwargs)


def split_into_datasets(data, split_key="split", **kwargs):
    datasets = {}

    if isinstance(data, dict):
        splits = list(data.keys())
    else:
        adata = data
        splits = list(adata.obs[split_key].unique())

    # ensure train comes first for encoder
    if "train" in splits:
        splits.remove("train")
        splits.insert(0, "train")

    label_encoder = None
    for i, split in enumerate(splits):
        if isinstance(data, dict):
            adata_split = data[split]
        else:
            adata_split = adata[
                adata.obs[split_key] == split
            ].copy()  # make copy to avoid slow lazy indexing

        datasets[split] = AnnDataSet(adata_split, label_encoder=label_encoder, **kwargs)

        if i == 0:
            label_encoder = datasets[split]._label_encoder

    return datasets
