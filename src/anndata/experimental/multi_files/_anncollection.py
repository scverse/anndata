from __future__ import annotations

import warnings
from collections.abc import Callable, Mapping
from functools import reduce
from itertools import chain, pairwise
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from h5py import Dataset

from ..._core.aligned_mapping import AxisArrays
from ..._core.anndata import AnnData
from ..._core.index import _normalize_index, _normalize_indices
from ..._core.merge import concat_arrays, inner_concat_aligned_mapping
from ..._core.sparse_dataset import BaseCompressedSparseDataset
from ..._core.views import _resolve_idx
from ...compat import _map_cat_to_str, old_positionals

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    from ..._core.index import Index

ATTRS = ["obs", "obsm", "layers"]


def _merge(arrs):
    rxers = [lambda x, fill_value, axis: x] * len(arrs)
    return concat_arrays(arrs, rxers)


def _select_convert(key, convert, arr=None):
    key_convert = None

    if callable(convert):
        key_convert = convert
    elif isinstance(convert, dict) and key in convert:
        key_convert = convert[key]

    if arr is not None:
        return key_convert(arr) if key_convert is not None else arr
    else:
        return key_convert


def _harmonize_types(attrs_keys, adatas):
    attrs_keys_types = {}

    def check_type(attr, key=None):
        arrs = []
        for a in adatas:
            attr_arr = getattr(a, attr)
            if key is not None:
                attr_arr = attr_arr[key]
            arrs.append(attr_arr)
        # hacky but numpy find_common_type doesn't work with categoricals
        try:
            dtype = _merge([arr[:1] for arr in arrs]).dtype
        except ValueError:
            dtype = _merge([arr[:1, :1] for arr in arrs]).dtype
        return dtype

    for attr, keys in attrs_keys.items():
        if len(keys) == 0:
            continue
        attrs_keys_types[attr] = {}
        for key in keys:
            attrs_keys_types[attr][key] = check_type(attr, key)

    attrs_keys_types["X"] = check_type("X")

    return attrs_keys_types


class _ConcatViewMixin:
    def _resolve_idx(self, oidx, vidx):
        adatas_oidx = []
        reverse = None

        old_oidx = getattr(self, "oidx", None)
        if old_oidx is not None:
            oidx = _resolve_idx(old_oidx, oidx, self.limits[-1])

        if isinstance(oidx, slice):
            start, stop, step = oidx.indices(self.limits[-1])
            oidx = np.arange(start, stop, step)
        else:
            oidx = np.array([oidx]) if isinstance(oidx, int) else oidx
        u_oidx = oidx

        if len(self.adatas) == 1:
            return [u_oidx], oidx, vidx, reverse

        iter_limits = list(pairwise(chain([0], self.limits)))

        n_adatas_used = 0
        for lower, upper in iter_limits:
            if np.any((u_oidx >= lower) & (u_oidx < upper)):
                n_adatas_used += 1

        need_reverse = (
            self.indices_strict
            and n_adatas_used > 1
            and u_oidx.size > 1
            and np.any(u_oidx[:-1] > u_oidx[1:])
        )
        if need_reverse:
            u_oidx, reverse = np.unique(u_oidx, return_inverse=True)

        for lower, upper in iter_limits:
            mask = (u_oidx >= lower) & (u_oidx < upper)
            adatas_oidx.append(u_oidx[mask] - lower if mask.any() else None)

        old_vidx = getattr(self, "vidx", None)
        if old_vidx is not None:
            vidx = _resolve_idx(old_vidx, vidx, self.adatas[0].n_vars)
        if isinstance(vidx, int):
            vidx = np.array([vidx])

        return adatas_oidx, oidx, vidx, reverse


class _IterateViewMixin:
    @old_positionals("axis", "shuffle", "drop_last")
    def iterate_axis(
        self,
        batch_size: int,
        *,
        axis: Literal[0, 1] = 0,
        shuffle: bool = False,
        drop_last: bool = False,
    ):
        """Iterate the lazy object over an axis.

        Parameters
        ----------
        batch_size
            How many samples to put into a batch when iterating.
        axis
            The axis to iterate over.
        shuffle
            Set to `True` to have the indices reshuffled before iterating.
        drop_last
            Set to `True` to drop a batch with the length lower than `batch_size`.
        """
        if axis not in {0, 1}:
            msg = "Axis should be either 0 or 1."
            raise ValueError(msg)

        n = self.shape[axis]
        indices = np.random.permutation(n).tolist() if shuffle else list(range(n))

        for i in range(0, n, batch_size):
            idx = indices[i : min(i + batch_size, n)]
            batch = self[:, idx] if axis == 1 else self[idx]
            # only happens if the last batch is smaller than batch_size
            if len(batch) < batch_size and drop_last:
                continue

            yield batch, idx


class MapObsView:
    def __init__(  # noqa: PLR0913
        self,
        attr,
        adatas,
        keys,
        *,
        adatas_oidx,
        adatas_vidx=None,
        convert=None,
        reverse=None,
        dtypes=None,
        obs_names=None,
    ):
        self.adatas = adatas
        self._keys = keys
        self.adatas_oidx = adatas_oidx
        self.adatas_vidx = adatas_vidx
        self.attr = attr
        self.convert = convert
        self.reverse = reverse
        self.dtypes = dtypes
        self.obs_names = obs_names

    def __getitem__(self, key: str, *, use_convert: bool = True):
        if self._keys is not None and key not in self._keys:
            msg = f"No {key} in {self.attr} view"
            raise KeyError(msg)

        arrs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is None:
                continue

            arr = getattr(self.adatas[i], self.attr)[key]
            vidx = self.adatas_vidx[i] if self.adatas_vidx is not None else None
            idx = (oidx, vidx) if vidx is not None else oidx

            if isinstance(arr, pd.DataFrame):
                arrs.append(arr.iloc[idx])
            else:
                if vidx is not None:
                    idx = np.ix_(*idx) if not isinstance(idx[1], slice) else idx
                arrs.append(arr.iloc[idx] if isinstance(arr, pd.Series) else arr[idx])

        if len(arrs) > 1:
            _arr = _merge(arrs)
            _arr = _arr if self.reverse is None else _arr[self.reverse]
        else:
            _arr = arrs[0]
            # what if it is a dataframe?
            if self.dtypes is not None:
                _arr = _arr.astype(self.dtypes[key], copy=False)

        if self.convert is not None and use_convert:
            _arr = _select_convert(key, self.convert, _arr)

        return _arr

    def keys(self):
        if self._keys is not None:
            return self._keys
        else:
            return list(getattr(self.adatas[0], self.attr).keys())

    @old_positionals("use_convert")
    def to_dict(self, keys: Iterable[str] | None = None, *, use_convert=True):
        dct = {}
        keys = self.keys() if keys is None else keys
        for key in keys:
            dct[key] = self.__getitem__(key, use_convert=use_convert)
        return dct

    @property
    def df(self):
        if self.attr != "obs":
            return None
        return pd.DataFrame(self.to_dict(use_convert=False), index=self.obs_names)

    def __repr__(self):
        descr = f"View of {self.attr} with keys: {str(self.keys())[1:-1]}"
        return descr


class AnnCollectionView(_ConcatViewMixin, _IterateViewMixin):
    """\
    An object to access the observation attributes of `adatas` in AnnCollection.

    Created as a result of subsetting an :class:`~anndata.experimental.AnnCollection` object.
    An object of this class can have `.obs`, `.obsm`, `.layers`, `.X` depending on the
    results of joins in the reference AnnCollection object.

    Notes
    -----
    Nothing is copied until keys of the attributes or `.X` are accessed.
    """

    def __init__(self, reference, convert, resolved_idx):
        self.reference = reference

        self.indices_strict = self.reference.indices_strict

        self.adatas = self.reference.adatas
        self.limits = self.reference.limits

        self.adatas_oidx, self.oidx, self.vidx, self.reverse = resolved_idx

        self.adatas_vidx = []

        for i, vidx in enumerate(self.reference.adatas_vidx):
            if vidx is None:
                self.adatas_vidx.append(self.vidx)
            else:
                new_vidx = _resolve_idx(vidx, self.vidx, self.adatas[i].n_vars)
                self.adatas_vidx.append(new_vidx)

        self._view_attrs_keys = self.reference._view_attrs_keys
        self._attrs = self.reference._attrs

        self._dtypes = self.reference._dtypes

        self._layers_view, self._obsm_view, self._obs_view = None, None, None
        self._X = None

        self._convert = None
        self._convert_X = None
        self.convert = convert

    def _lazy_init_attr(self, attr: str, *, set_vidx: bool = False):
        if getattr(self, f"_{attr}_view") is not None:
            return
        keys = None
        attr_dtypes = None
        if attr in self._view_attrs_keys:
            reverse = self.reverse
            keys = self._view_attrs_keys[attr]
            if len(keys) == 0:
                return
            adatas = self.adatas
            adatas_oidx = self.adatas_oidx
            if self._dtypes is not None:
                attr_dtypes = self._dtypes[attr]
        else:
            reverse = None
            adatas = [self.reference]
            adatas_oidx = [self.oidx]
        adatas_vidx = self.adatas_vidx if set_vidx else None

        attr_convert = None
        if self.convert is not None:
            attr_convert = _select_convert(attr, self.convert)

        obs_names = self.obs_names if attr == "obs" else None

        setattr(
            self,
            f"_{attr}_view",
            MapObsView(
                attr,
                adatas,
                keys,
                adatas_oidx=adatas_oidx,
                adatas_vidx=adatas_vidx,
                convert=attr_convert,
                reverse=reverse,
                dtypes=attr_dtypes,
                obs_names=obs_names,
            ),
        )

    def _gather_X(self):
        if self._X is not None:
            return self._X

        Xs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is None:
                continue

            adata = self.adatas[i]
            X = adata.X
            vidx = self.adatas_vidx[i]

            if isinstance(X, Dataset):
                reverse = None
                if oidx.size > 1 and np.any(oidx[:-1] >= oidx[1:]):
                    oidx, reverse = np.unique(oidx, return_inverse=True)

                # TODO: fix memory inefficient approach of X[oidx][:, vidx]
                arr = X[oidx, vidx] if isinstance(vidx, slice) else X[oidx][:, vidx]
                Xs.append(arr if reverse is None else arr[reverse])
            elif isinstance(X, BaseCompressedSparseDataset):
                # very slow indexing with two arrays
                if isinstance(vidx, slice) or len(vidx) <= 1000:
                    Xs.append(X[oidx, vidx])
                else:
                    Xs.append(X[oidx][:, vidx])
            else:
                # if vidx is present it is less memory efficient
                idx = oidx, vidx
                idx = np.ix_(*idx) if not isinstance(vidx, slice) else idx
                Xs.append(X[idx])

        if len(Xs) > 1:
            _X = _merge(Xs)
            # todo: get rid of reverse for dense arrays
            _X = _X if self.reverse is None else _X[self.reverse]
        else:
            _X = Xs[0]
            if self._dtypes is not None:
                _X = _X.astype(self._dtypes["X"], copy=False)

        self._X = _X

        return _X

    @property
    def X(self):
        """Lazy subset of data matrix.

        The data matrix formed from the `.X` attributes of the underlying `adatas`,
        properly reindexed and lazily merged.
        Nothing is copied until `.X` is accessed, no real concatenation of the
        underlying `.X` attributes is done.
        """
        # inconsistent behavior here, _X can be changed,
        # but the other attributes can't be changed.
        # maybe do return ... _X.copy() or _X.setflags(write=False)

        _X = self._gather_X()

        return self._convert_X(_X) if self._convert_X is not None else _X

    @property
    def layers(self):
        """Lazy subset of layers.

        The layers attribute formed from lazy inner join and subsetting of the `.layers`
        of the underlying `adatas`. No copy is made until you access a key from `.layers`,
        only the subset of the accessed key is copied.

        To get `.layers` as a dictionary, use `.layers.to_dict()`. You can also specify keys
        to include in the dict `.layers.to_dict(keys=['key1', 'key2'])` and if you want
        converters to be turned off when copying to dict `.layers.to_dict(use_convert=False)`.
        """
        self._lazy_init_attr("layers", set_vidx=True)
        return self._layers_view

    @property
    def obsm(self):
        """Lazy subset of multi-dimensional annotation of observations.

        Points to the `.obsm` attributes of the underlying adatas to `.obsm` of the parent
        AnnCollection object depending on the `join_obsm` option of the AnnCollection object.
        See the docs of :class:`~anndata.experimental.AnnCollection` for details.
        Copy rules are the same as for `.layers`, i.e. everything is lazy.

        To get `.obsm` as a dictionary, use `.obsm.to_dict()`. You can also specify keys
        to include in the dict `.obsm.to_dict(keys=['key1', 'key2'])` and if you want
        converters to be turned off when copying to dict `.obsm.to_dict(use_convert=False)`.
        """
        self._lazy_init_attr("obsm")
        return self._obsm_view

    @property
    def obs(self):
        """Lazy suset of one-dimensional annotation of observations.

        Points to the `.obs` attributes of the underlying adatas to `.obs` of the parent
        AnnCollection object depending on the `join_obs` option of the AnnCollection object.
        See the docs of `~anndata.experimental.AnnCollection` for details.
        Copy rules are the same as for `.layers`, i.e. everything is lazy.

        To get `.obs` as a DataFrame, use `.obs.df`.
        To get `.obs` as a dictionary, use `.obs.to_dict()`. You can also specify keys
        to include in the dict `.obs.to_dict(keys=['key1', 'key2'])` and if you want
        converters to be turned off when copying to dict `.obs.to_dict(use_convert=False)`.
        """
        self._lazy_init_attr("obs")
        return self._obs_view

    @property
    def obs_names(self):
        """Names of observations of this subset object."""
        return self.reference.obs_names[self.oidx]

    @property
    def var_names(self):
        """Names of variables of this subset object."""
        return self.reference.var_names[self.vidx]

    @property
    def shape(self):
        """Shape of the lazily concatenated subset of the data matrix."""
        return len(self.obs_names), len(self.var_names)

    @property
    def n_obs(self):
        """Number of observations."""
        return self.shape[0]

    @property
    def n_vars(self):
        """Number of variables/features."""
        return self.shape[1]

    @property
    def convert(self):
        """On the fly converters for keys of attributes and data matrix.

        A function or a Mapping of functions which will be applied
        to the values of attributes (`.X`) or to specific keys of these attributes
        (`.obs`, `.obsm`, `.layers`).
        The keys of the Mapping should correspond to the attributes or keys of the
        attributes (hierarchically) and the values should be functions used for conversion.

        Examples
        ----------
        ::

            {
                # densify .X
                "X": lambda a: a.toarray() if issparse(a) else a,
                # change dtype for all keys of .obsm
                "obsm": lambda a: np.asarray(a, dtype="float32"),
                # change type only for one key of .obs
                "obs": dict(key1=lambda c: c.astype(str)),
            }
        """
        return self._convert

    @convert.setter
    def convert(self, value):
        self._convert = value
        self._convert_X = _select_convert("X", self._convert)
        for attr in ATTRS:
            setattr(self, f"_{attr}_view", None)

    def __len__(self):
        return len(self.obs_names)

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnCollectionView(self.reference, self.convert, resolved_idx)

    @property
    def has_backed(self):
        """`True` if the current subset of `adatas` has backed objects, `False` otherwise."""
        for i, adata in enumerate(self.adatas):
            if adata.isbacked and self.adatas_oidx[i] is not None:
                return True
        return False

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnCollectionView object with n_obs × n_vars = {n_obs} × {n_vars}"
        all_attrs_keys = self._view_attrs_keys.copy()
        for attr in self._attrs:
            all_attrs_keys[attr] = list(getattr(self.reference, attr).keys())
        for attr, keys in all_attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(keys)[1:-1]}"
        return descr

    @old_positionals("ignore_X", "ignore_layers")
    def to_adata(self, *, ignore_X: bool = False, ignore_layers: bool = False):
        """Convert this AnnCollectionView object to an AnnData object.

        Parameters
        ----------
        ignore_X
            if `True`, adds `.X` to the AnnData object.
        ignore_layers
            if `True`, copies `.layers` to the AnnData object.
        """
        if ignore_layers or self.layers is None:
            layers = None
        else:
            layers = self.layers.to_dict(use_convert=False)
        obsm = None if self.obsm is None else self.obsm.to_dict(use_convert=False)
        obs = (
            None
            if self.obs is None
            else pd.DataFrame(self.obs.to_dict(use_convert=False))
        )

        if ignore_X:
            X = None
            shape = self.shape
        else:
            X = self._gather_X()
            shape = None

        adata = AnnData(X, obs=obs, obsm=obsm, layers=layers, shape=shape)
        adata.obs_names = self.obs_names
        adata.var_names = self.var_names
        return adata

    @property
    def attrs_keys(self):
        """Dict of all accessible attributes and their keys."""
        return self.reference.attrs_keys


DictCallable = dict[str, Callable]
ConvertType = Callable | dict[str, Callable | DictCallable]


class AnnCollection(_ConcatViewMixin, _IterateViewMixin):
    """\
    Lazily concatenate AnnData objects along the `obs` axis.

    This class doesn't copy data from underlying AnnData objects, but lazily subsets using a joint
    index of observations and variables. It also allows on-the-fly application of prespecified
    converters to `.obs` attributes of the AnnData objects.

    Subsetting of this object returns an `AnnCollectionView`, which provides views of `.obs`,
    `.obsm`, `.layers`, `.X` from the underlying AnnData objects.

    Parameters
    ----------
    adatas
        The objects to be lazily concatenated.
        If a Mapping is passed, keys are used for the `keys` argument and values are concatenated.
    join_obs
        If "inner" specified all `.obs` attributes from `adatas` will be inner joined
        and copied to this object.
        If "outer" specified all `.obsm` attributes from `adatas` will be outer joined
        and copied to this object.
        For "inner" and "outer" subset objects will access `.obs` of this object,
        not the original `.obs` attributes of `adatas`.
        If `None`, nothing is copied to this object's `.obs`, a subset object will directly
        access `.obs` attributes of `adatas` (with proper reindexing and dtype conversions).
        For `None`the inner join rule is used to select columns of `.obs` of `adatas`.
    join_obsm
        If "inner" specified all `.obsm` attributes from `adatas` will be inner joined
        and copied to this object. Subset objects will access `.obsm` of this object,
        not the original `.obsm` attributes of `adatas`.
        If `None`, nothing is copied to this object's `.obsm`, a subset object will directly
        access `.obsm` attributes of `adatas` (with proper reindexing and dtype conversions).
        For both options the inner join rule for the underlying `.obsm` attributes is used.
    join_vars
        Specify how to join `adatas` along the var axis. If `None`, assumes all `adatas`
        have the same variables. If "inner", the intersection of all variables in
        `adatas` will be used.
    label
        Column in `.obs` to place batch information in.
        If it's None, no column is added.
    keys
        Names for each object being added. These values are used for column values for
        `label` or appended to the index if `index_unique` is not `None`. Defaults to
        incrementing integer labels.
    index_unique
        Whether to make the index unique by using the keys. If provided, this
        is the delimiter between "{orig_idx}{index_unique}{key}". When `None`,
        the original indices are kept.
    convert
        You can pass a function or a Mapping of functions which will be applied
        to the values of attributes (`.obs`, `.obsm`, `.layers`, `.X`) or to specific
        keys of these attributes in the subset object.
        Specify an attribute and a key (if needed) as keys of the passed Mapping
        and a function to be applied as a value.
    harmonize_dtypes
        If `True`, all retrieved arrays from subset objects will have the same dtype.
    indices_strict
        If  `True`, arrays from the subset objects will always have the same order
        of indices as in selection used to subset.
        This parameter can be set to `False` if the order in the returned arrays
        is not important, for example, when using them for stochastic gradient descent.
        In this case the performance of subsetting can be a bit better.

    Examples
    ----------
    >>> from scanpy.datasets import pbmc68k_reduced, pbmc3k_processed
    >>> adata1, adata2 = pbmc68k_reduced(), pbmc3k_processed()
    >>> adata1.shape
    (700, 765)
    >>> adata2.shape
    (2638, 1838)
    >>> dc = AnnCollection([adata1, adata2], join_vars='inner')
    >>> dc
    AnnCollection object with n_obs × n_vars = 3338 × 208
      constructed from 2 AnnData objects
        view of obsm: 'X_pca', 'X_umap'
        obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain'
    >>> batch = dc[100:200] # AnnCollectionView
    >>> batch
    AnnCollectionView object with n_obs × n_vars = 100 × 208
        obsm: 'X_pca', 'X_umap'
        obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain'
    >>> batch.X.shape
    (100, 208)
    >>> len(batch.obs['louvain'])
    100
    """

    @old_positionals(
        "join_obs",
        "join_obsm",
        "join_vars",
        "label",
        "keys",
        "index_unique",
        "convert",
        "harmonize_dtypes",
        "indices_strict",
    )
    def __init__(  # noqa: PLR0912, PLR0913, PLR0915
        self,
        adatas: Sequence[AnnData] | dict[str, AnnData],
        *,
        join_obs: Literal["inner", "outer"] | None = "inner",
        join_obsm: Literal["inner"] | None = None,
        join_vars: Literal["inner"] | None = None,
        label: str | None = None,
        keys: Sequence[str] | None = None,
        index_unique: str | None = None,
        convert: ConvertType | None = None,
        harmonize_dtypes: bool = True,
        indices_strict: bool = True,
    ):
        if isinstance(adatas, Mapping):
            if keys is not None:
                msg = (
                    "Cannot specify categories in both mapping keys and using `keys`. "
                    "Only specify this once."
                )
                raise TypeError(msg)
            keys, adatas = list(adatas.keys()), list(adatas.values())
        else:
            adatas = list(adatas)

        # check if the variables are the same in all adatas
        self.adatas_vidx = [None for adata in adatas]
        vars_names_list = [adata.var_names for adata in adatas]
        vars_eq = all(adatas[0].var_names.equals(vrs) for vrs in vars_names_list[1:])
        if vars_eq:
            self.var_names = adatas[0].var_names
        elif join_vars == "inner":
            var_names = reduce(pd.Index.intersection, vars_names_list)
            self.adatas_vidx = []
            for adata in adatas:
                if var_names.equals(adata.var_names):
                    self.adatas_vidx.append(None)
                else:
                    adata_vidx = _normalize_index(var_names, adata.var_names)
                    self.adatas_vidx.append(adata_vidx)
            self.var_names = var_names
        else:
            msg = (
                "Adatas have different variables. "
                "Please specify join_vars='inner' for intersection."
            )
            raise ValueError(msg)

        concat_indices = pd.concat(
            [pd.Series(a.obs_names) for a in adatas], ignore_index=True
        )
        if keys is None:
            keys = np.arange(len(adatas)).astype(str)
        label_col = pd.Categorical.from_codes(
            np.repeat(np.arange(len(adatas)), [a.shape[0] for a in adatas]),
            categories=keys,
        )
        if index_unique is not None:
            concat_indices = concat_indices.str.cat(
                _map_cat_to_str(label_col), sep=index_unique
            )
        self.obs_names = pd.Index(concat_indices)

        if not self.obs_names.is_unique:
            msg = "Observation names are not unique."
            warnings.warn(msg, UserWarning, stacklevel=2)

        view_attrs = ATTRS.copy()

        self._attrs = []
        # process obs joins
        if join_obs is not None:
            view_attrs.remove("obs")
            self._attrs.append("obs")
            concat_annot = pd.concat(
                [a.obs for a in adatas], join=join_obs, ignore_index=True
            )
            concat_annot.index = self.obs_names
            self._obs = concat_annot
        else:
            self._obs = pd.DataFrame(index=self.obs_names)

        if label is not None:
            self._obs[label] = label_col

        # process obsm inner join
        self._obsm = None
        if join_obsm == "inner":
            view_attrs.remove("obsm")
            self._attrs.append("obsm")
            self._obsm = inner_concat_aligned_mapping(
                [a.obsm for a in adatas], index=self.obs_names
            )
            self._obsm = (
                AxisArrays(self, axis=0, store={}) if self._obsm == {} else self._obsm
            )

        # process inner join of views
        self._view_attrs_keys = {}
        for attr in view_attrs:
            self._view_attrs_keys[attr] = list(getattr(adatas[0], attr).keys())

        for a in adatas[1:]:
            for attr, keys in self._view_attrs_keys.items():
                ai_attr = getattr(a, attr)
                a0_attr = getattr(adatas[0], attr)
                new_keys = []
                for key in keys:
                    if key in ai_attr:
                        a0_ashape = a0_attr[key].shape
                        ai_ashape = ai_attr[key].shape
                        if (
                            len(a0_ashape) < 2
                            or a0_ashape[1] == ai_ashape[1]
                            or attr == "layers"
                        ):
                            new_keys.append(key)
                self._view_attrs_keys[attr] = new_keys

        self.adatas = adatas

        self.limits = [adatas[0].n_obs]
        for i in range(len(adatas) - 1):
            self.limits.append(self.limits[i] + adatas[i + 1].n_obs)

        # init converter
        self._convert = convert

        self._dtypes = None
        if len(adatas) > 1 and harmonize_dtypes:
            self._dtypes = _harmonize_types(self._view_attrs_keys, self.adatas)

        self.indices_strict = indices_strict

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnCollectionView(self, self.convert, resolved_idx)

    @property
    def convert(self):
        """On the fly converters for keys of attributes and data matrix.

        A function or a Mapping of functions which will be applied
        to the values of attributes (`.X`) or to specific keys of these attributes
        (`.obs`, `.obsm`, `.layers`) of subset objects. The converters are not
        applied to `.obs` and `.obsm` (if present) of this object, only to the attributes
        of subset objects.
        The keys of the Mapping should correspond to the attributes or keys of the
        attributes (hierarchically) and the values should be functions used for conversion.

        Examples
        --------
        ::

            {
                # densify .X
                "X": lambda a: a.toarray() if issparse(a) else a,
                # change dtype for all keys of .obsm
                "obsm": lambda a: np.asarray(a, dtype="float32"),
                # change type only for one key of .obs
                "obs": dict(key1=lambda c: c.astype(str)),
            }
        """
        return self._convert

    @convert.setter
    def convert(self, value):
        self._convert = value

    @property
    def obs(self):
        """One-dimensional annotation of observations.

        If `join_obs` was set to "inner" and "outer", subset objects' `.obs`
        will point to this `.obs`; otherwise, to `.obs` of the underlying objects (`adatas`).
        """
        return self._obs

    @property
    def obsm(self):
        """Multi-dimensional annotation of observations.

        If `join_obsm` was set to "inner", subset objects' `.obsm`
        will point to this `.obsm`; otherwise, to `.obsm` of the underlying objects (`adatas`).
        In the latter case, `.obsm` of this object will be `None`.
        """
        return self._obsm

    @property
    def shape(self):
        """Shape of the lazily concatenated data matrix"""
        return self.limits[-1], len(self.var_names)

    @property
    def n_obs(self):
        """Number of observations."""
        return self.shape[0]

    @property
    def n_vars(self):
        """Number of variables/features."""
        return self.shape[1]

    def __len__(self):
        return self.limits[-1]

    def to_adata(self):
        """Convert this AnnCollection object to an AnnData object.

        The AnnData object won't have `.X`, only `.obs` and `.obsm`.
        """
        if "obs" in self._view_attrs_keys or "obsm" in self._view_attrs_keys:
            concat_view = self[self.obs_names]

        if "obsm" in self._view_attrs_keys:
            obsm = (
                concat_view.obsm.to_dict(use_convert=False)
                if concat_view.obsm is not None
                else None
            )
        else:
            obsm = self.obsm.copy()

        obs = self.obs.copy()
        if "obs" in self._view_attrs_keys and concat_view.obs is not None:
            for key, value in concat_view.obs.to_dict(use_convert=False).items():
                obs[key] = value

        adata = AnnData(X=None, obs=obs, obsm=obsm, shape=self.shape)
        adata.obs_names = self.obs_names
        adata.var_names = self.var_names
        return adata

    def lazy_attr(self, attr, key=None):
        """Get a subsettable key from an attribute (array-like) or an attribute.

        Returns a LazyAttrData object which provides subsetting over the specified
        attribute (`.obs` or `.obsm`) or over a key from this attribute.
        In the latter case, it acts as a lazy array.
        """
        return LazyAttrData(self, attr, key)

    @property
    def has_backed(self):
        """`True` if `adatas` have backed AnnData objects, `False` otherwise."""
        return any(adata.isbacked for adata in self.adatas)

    @property
    def attrs_keys(self):
        """Dict of all accessible attributes and their keys."""
        _attrs_keys = {}
        for attr in self._attrs:
            keys = list(getattr(self, attr).keys())
            _attrs_keys[attr] = keys
        _attrs_keys.update(self._view_attrs_keys)
        return _attrs_keys

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnCollection object with n_obs × n_vars = {n_obs} × {n_vars}"
        descr += f"\n  constructed from {len(self.adatas)} AnnData objects"
        for attr, keys in self._view_attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    view of {attr}: {str(keys)[1:-1]}"
        for attr in self._attrs:
            keys = list(getattr(self, attr).keys())
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(keys)[1:-1]}"
        if "obs" in self._view_attrs_keys:
            keys = list(self.obs.keys())
            if len(keys) > 0:
                descr += f"\n    own obs: {str(keys)[1:-1]}"

        return descr


class LazyAttrData(_IterateViewMixin):
    def __init__(self, adset: AnnCollection, attr: str, key: str | None = None):
        self.adset = adset
        self.attr = attr
        self.key = key

    def __getitem__(self, index):
        oidx = None
        vidx = None

        if isinstance(index, tuple) and self.attr in {"obs", "obsm"}:
            oidx = index[0]
            if len(index) > 1:
                vidx = index[1]

        view = self.adset[index] if oidx is None else self.adset[oidx]
        attr_arr = getattr(view, self.attr)
        if self.key is not None:
            attr_arr = attr_arr[self.key]
        return attr_arr if vidx is None else attr_arr[:, vidx]

    @property
    def shape(self):
        shape = self.adset.shape
        if self.attr in {"X", "layers"}:
            return shape
        elif self.attr == "obs":
            return (shape[0],)
        elif self.attr == "obsm" and self.key is not None:
            return shape[0], self[:1].shape[1]
        else:
            return None

    @property
    def ndim(self):
        return len(self.shape) if self.shape is not None else 0

    @property
    def dtype(self):
        _dtypes = self.adset._dtypes
        if _dtypes is not None and self.attr in _dtypes:
            return _dtypes[self.attr][self.key]

        attr = self[:1]
        if hasattr(attr, "dtype"):
            return attr.dtype
        else:
            return None
