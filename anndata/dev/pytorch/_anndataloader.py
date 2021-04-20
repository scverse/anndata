from scipy.sparse import issparse
from ..._core.anndata import AnnData
from ..multi_files._anndataset import AnnDataConcatObs, _ConcatViewMixin
import numpy as np
import warnings

try:
    import torch
    from torch.utils.data import Sampler, Dataset, DataLoader
except ImportError:
    warnings.warn("Сould not load pytorch.")
    Sampler, Dataset, DataLoader = object, object, object


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


def default_converter(arr, use_cuda):
    if np.issubdtype(arr.dtype, np.number):
        if issparse(arr):
            arr = arr.toarray()
        return torch.cuda.FloatTensor(arr) if use_cuda else torch.FloatTensor(arr)
    else:
        return arr


# AnnDataLoader has the same arguments as DataLoader, but uses BatchIndexSampler by default
class AnnDataLoader(DataLoader):
    def __init__(self, adata, batch_size=1, shuffle=False, **kwargs):

        if isinstance(adata, AnnData):
            join_obs = kwargs.pop("join_obs", "inner")
            join_obsm = kwargs.pop("join_obsm", None)
            label = kwargs.pop("label", None)
            keys = kwargs.pop("keys", None)
            index_unique = kwargs.pop("index_unique", None)

            use_cuda = kwargs.pop("use_cuda", False)

            _converter = lambda arr: default_converter(arr, use_cuda)
            convert = kwargs.pop("convert", _converter)

            dataset = AnnDataConcatObs(
                [adata],
                join_obs=join_obs,
                join_obsm=join_obsm,
                label=label,
                keys=keys,
                index_unique=index_unique,
                convert=convert,
            )

        elif isinstance(adata, _ConcatViewMixin):
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
