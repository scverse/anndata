from __future__ import annotations

from collections.abc import Mapping
from copy import copy
from functools import partial
from importlib.util import find_spec
from math import ceil
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import issparse

from ..._core.anndata import AnnData
from ...compat import old_positionals
from ..multi_files._anncollection import AnnCollection, _ConcatViewMixin

if find_spec("torch") or TYPE_CHECKING:
    import torch
    from torch.utils.data import BatchSampler, DataLoader, Sampler
else:
    Sampler, BatchSampler, DataLoader = object, object, object

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Sequence
    from typing import TypeAlias, Union

    from scipy.sparse import spmatrix

    # need to use Union because of autodoc_mock_imports
    Array: TypeAlias = Union[torch.Tensor, np.ndarray, spmatrix]  # noqa: UP007


# Custom sampler to get proper batches instead of joined separate indices
# maybe move to multi_files
class BatchIndexSampler(Sampler):
    @old_positionals("batch_size", "shuffle", "drop_last")
    def __init__(
        self,
        n_obs: int,
        *,
        batch_size: int,
        shuffle: bool = False,
        drop_last: bool = False,
    ) -> None:
        self.n_obs = n_obs
        self.batch_size = batch_size if batch_size < n_obs else n_obs
        self.shuffle = shuffle
        self.drop_last = drop_last

    def __iter__(self) -> Generator[list[int], None, None]:
        indices: list[int]
        if self.shuffle:
            indices = np.random.permutation(self.n_obs).tolist()
        else:
            indices = list(range(self.n_obs))

        for i in range(0, self.n_obs, self.batch_size):
            batch = indices[i : min(i + self.batch_size, self.n_obs)]

            # only happens if the last batch is smaller than batch_size
            if len(batch) < self.batch_size and self.drop_last:
                continue

            yield batch

    def __len__(self) -> int:
        if self.drop_last:
            length = self.n_obs // self.batch_size
        else:
            length = ceil(self.n_obs / self.batch_size)

        return length


# maybe replace use_cuda with explicit device option
def default_converter(arr: Array, *, use_cuda: bool, pin_memory: bool):
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


def _convert_on_top(
    convert: Callable[[Array], Array] | None | Mapping[str, Callable[[Array], Array]],
    top_convert: Callable[[Array], Array],
    attrs_keys: Sequence[str] | Mapping[str, Sequence[str]],
):
    if convert is None:
        new_convert = top_convert
    elif callable(convert):

        def compose_convert(arr):
            return top_convert(convert(arr))

        new_convert = compose_convert
    else:
        new_convert = {}
        for attr in attrs_keys:
            if attr not in convert:
                new_convert[attr] = top_convert
            else:
                as_ks: Sequence[str] | None
                if not isinstance(attrs_keys, Mapping):
                    as_ks = None
                else:
                    as_ks = attrs_keys[attr]
                new_convert[attr] = _convert_on_top(convert[attr], top_convert, as_ks)
    return new_convert


# AnnLoader has the same arguments as DataLoader, but uses BatchIndexSampler by default
class AnnLoader(DataLoader):
    """\
    PyTorch DataLoader for AnnData objects.

    Builds DataLoader from a sequence of AnnData objects, from an
    :class:`~anndata.experimental.AnnCollection` object or from an `AnnCollectionView` object.
    Takes care of the required conversions.

    Parameters
    ----------
    adatas
        `AnnData` objects or an `AnnCollection` object from which to load the data.
    batch_size
        How many samples per batch to load.
    shuffle
        Set to `True` to have the data reshuffled at every epoch.
    use_default_converter
        Use the default converter to convert arrays to pytorch tensors, transfer to
        the default cuda device (if `use_cuda=True`), do memory pinning (if `pin_memory=True`).
        If you pass an AnnCollection object with prespecified converters, the default converter
        won't overwrite these converters but will be applied on top of them.
    use_cuda
        Transfer pytorch tensors to the default cuda device after conversion.
        Only works if `use_default_converter=True`
    **kwargs
        Arguments for PyTorch DataLoader. If `adatas` is not an `AnnCollection` object, then also
        arguments for `AnnCollection` initialization.
    """

    @old_positionals("batch_size", "shuffle", "use_default_converter", "use_cuda")
    def __init__(
        self,
        adatas: Sequence[AnnData] | dict[str, AnnData],
        *,
        batch_size: int = 1,
        shuffle: bool = False,
        use_default_converter: bool = True,
        use_cuda: bool = False,
        **kwargs,
    ):
        if isinstance(adatas, AnnData):
            adatas = [adatas]

        if isinstance(adatas, list | tuple | dict):
            join_obs = kwargs.pop("join_obs", "inner")
            join_obsm = kwargs.pop("join_obsm", None)
            label = kwargs.pop("label", None)
            keys = kwargs.pop("keys", None)
            index_unique = kwargs.pop("index_unique", None)
            convert = kwargs.pop("convert", None)
            harmonize_dtypes = kwargs.pop("harmonize_dtypes", True)
            indices_strict = kwargs.pop("indices_strict", True)

            dataset = AnnCollection(
                adatas,
                join_obs=join_obs,
                join_obsm=join_obsm,
                label=label,
                keys=keys,
                index_unique=index_unique,
                convert=convert,
                harmonize_dtypes=harmonize_dtypes,
                indices_strict=indices_strict,
            )

        elif isinstance(adatas, _ConcatViewMixin):
            dataset = copy(adatas)
        else:
            msg = "adata should be of type AnnData or AnnCollection."
            raise ValueError(msg)

        if use_default_converter:
            pin_memory = kwargs.pop("pin_memory", False)
            _converter = partial(
                default_converter, use_cuda=use_cuda, pin_memory=pin_memory
            )
            dataset.convert = _convert_on_top(
                dataset.convert, _converter, dict(dataset.attrs_keys, X=[])
            )

        has_sampler = "sampler" in kwargs
        has_batch_sampler = "batch_sampler" in kwargs

        has_worker_init_fn = (
            "worker_init_fn" in kwargs and kwargs["worker_init_fn"] is not None
        )
        has_workers = "num_workers" in kwargs and kwargs["num_workers"] > 0
        use_parallel = has_worker_init_fn or has_workers

        if (
            batch_size is not None
            and batch_size > 1
            and not has_batch_sampler
            and not use_parallel
        ):
            drop_last = kwargs.pop("drop_last", False)

            if has_sampler:
                sampler = kwargs.pop("sampler")
                sampler = BatchSampler(
                    sampler, batch_size=batch_size, drop_last=drop_last
                )
            else:
                sampler = BatchIndexSampler(
                    len(dataset),
                    batch_size=batch_size,
                    shuffle=shuffle,
                    drop_last=drop_last,
                )

            super().__init__(dataset, batch_size=None, sampler=sampler, **kwargs)
        else:
            super().__init__(dataset, batch_size=batch_size, shuffle=shuffle, **kwargs)
