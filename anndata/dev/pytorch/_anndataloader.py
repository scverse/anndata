from scipy.sparse import issparse
from ..._core.anndata import AnnData
from ..multi_files._anndataset import AnnDataSet, _ConcatViewMixin
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


# maybe replace use_cuda with explicit device option
def default_converter(arr, use_cuda, pin_memory):
    if isinstance(arr, torch.Tensor):
        if use_cuda:
            arr = arr.cuda()
        elif pin_memory:
            arr = arr.pin_memory()
    elif arr.dtype.name != "category" and np.issubdtype(arr.dtype, np.number):
        if issparse(arr):
            arr = arr.toarray()
        if use_cuda:
            arr = torch.tensor(arr, device="cuda")
        else:
            arr = torch.tensor(arr)
            arr = arr.pin_memory() if pin_memory else arr
    return arr


def _convert_on_top(convert, top_convert, attrs_keys):
    if convert is None:
        new_convert = top_convert
    elif callable(convert):
        new_convert = lambda arr: top_convert(convert(arr))
    else:
        new_convert = {}
        for attr in attrs_keys:
            if attr not in convert:
                new_convert[attr] = top_convert
            else:
                if isinstance(attrs_keys, list):
                    as_ks = None
                else:
                    as_ks = attrs_keys[attr]
                new_convert[attr] = _convert_on_top(convert[attr], top_convert, as_ks)
    return new_convert


# AnnDataLoader has the same arguments as DataLoader, but uses BatchIndexSampler by default
class AnnDataLoader(DataLoader):
    def __init__(
        self,
        adatas,
        batch_size=1,
        shuffle=False,
        use_default_converter=True,
        use_cuda=False,
        **kwargs,
    ):

        if isinstance(adatas, AnnData):
            adatas = [adatas]

        if isinstance(adatas, list) or isinstance(adatas, tuple):
            join_obs = kwargs.pop("join_obs", "inner")
            join_obsm = kwargs.pop("join_obsm", None)
            label = kwargs.pop("label", None)
            keys = kwargs.pop("keys", None)
            index_unique = kwargs.pop("index_unique", None)
            convert = kwargs.pop("convert", None)

            dataset = AnnDataSet(
                adatas,
                join_obs=join_obs,
                join_obsm=join_obsm,
                label=label,
                keys=keys,
                index_unique=index_unique,
                convert=convert,
            )

        elif isinstance(adatas, _ConcatViewMixin):
            dataset = adatas
        else:
            raise ValueError("adata should be of type AnnData or AnnDataSet.")

        if use_default_converter:
            pin_memory = kwargs.pop("pin_memory", False)
            _converter = lambda arr: default_converter(arr, use_cuda, pin_memory)
            dataset.convert = _convert_on_top(
                dataset.convert, _converter, dict(dataset.attrs_keys, X=[])
            )

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
