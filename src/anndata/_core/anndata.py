"""\
Main class and helper functions.
"""

from __future__ import annotations

import warnings
from collections import OrderedDict
from collections.abc import Mapping, MutableMapping, Sequence
from copy import copy, deepcopy
from functools import partial, singledispatchmethod
from pathlib import Path
from textwrap import dedent
from typing import TYPE_CHECKING, cast

import h5py
import numpy as np
import pandas as pd
from natsort import natsorted
from numpy import ma
from pandas.api.types import infer_dtype
from scipy import sparse
from scipy.sparse import issparse

from anndata._warnings import ImplicitModificationWarning

from .. import utils
from .._settings import settings
from ..compat import CSArray, DaskArray, ZarrArray, _move_adj_mtx, old_positionals
from ..logging import anndata_logger as logger
from ..utils import (
    axis_len,
    deprecated,
    ensure_df_homogeneous,
    raise_value_error_if_multiindex_columns,
)
from .access import ElementRef
from .aligned_df import _gen_dataframe
from .aligned_mapping import AlignedMappingProperty, AxisArrays, Layers, PairwiseArrays
from .file_backing import AnnDataFileManager, to_memory
from .index import _normalize_indices, _subset, get_vector
from .raw import Raw
from .sparse_dataset import BaseCompressedSparseDataset, sparse_dataset
from .storage import coerce_array
from .views import (
    DictView,
    _resolve_idxs,
    as_view,
)
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Iterable
    from os import PathLike
    from typing import Any, ClassVar, Literal

    from zarr.storage import StoreLike

    from ..compat import Index1D, XDataset
    from ..typing import XDataType
    from .aligned_mapping import AxisArraysView, LayersView, PairwiseArraysView
    from .index import Index


class AnnData(metaclass=utils.DeprecationMixinMeta):  # noqa: PLW1641
    """\
    An annotated data matrix.

    .. figure:: ../_static/img/anndata_schema.svg
       :width: 260px
       :align: right
       :class: dark-light

    :class:`~anndata.AnnData` stores a data matrix :attr:`X` together with annotations
    of observations :attr:`obs` (:attr:`obsm`, :attr:`obsp`),
    variables :attr:`var` (:attr:`varm`, :attr:`varp`),
    and unstructured annotations :attr:`uns`.

    An :class:`~anndata.AnnData` object `adata` can be sliced like a
    :class:`~pandas.DataFrame`,
    for instance `adata_subset = adata[:, list_of_variable_names]`.
    :class:`~anndata.AnnData`’s basic structure is similar to R’s ExpressionSet
    [Huber15]_. If setting an `.h5ad`-formatted HDF5 backing file `.filename`,
    data remains on the disk but is automatically loaded into memory if needed.

    Parameters
    ----------
    X
        A #observations × #variables data matrix. A view of the data is used if the
        data type matches, otherwise, a copy is made.
    obs
        Key-indexed one-dimensional observations annotation of length #observations.
    var
        Key-indexed one-dimensional variables annotation of length #variables.
    uns
        Key-indexed unstructured annotation.
    obsm
        Key-indexed multi-dimensional observations annotation of length #observations.
        If passing a :class:`~numpy.ndarray`, it needs to have a structured datatype.
    varm
        Key-indexed multi-dimensional variables annotation of length #variables.
        If passing a :class:`~numpy.ndarray`, it needs to have a structured datatype.
    layers
        Key-indexed multi-dimensional arrays aligned to dimensions of `X`.
    shape
        Shape tuple (#observations, #variables). Can only be provided if `X` is `None`.
    filename
        Name of backing file. See :class:`h5py.File`.
    filemode
        Open mode of backing file. See :class:`h5py.File`.

    See Also
    --------
    io.read_h5ad
    io.read_csv
    io.read_excel
    io.read_hdf
    io.read_loom
    io.read_zarr
    io.read_mtx
    io.read_text
    io.read_umi_tools

    Notes
    -----
    :class:`~anndata.AnnData` stores observations (samples) of variables/features
    in the rows of a matrix.
    This is the convention of the modern classics of statistics [Hastie09]_
    and machine learning [Murphy12]_,
    the convention of dataframes both in R and Python and the established statistics
    and machine learning packages in Python (statsmodels_, scikit-learn_).

    Single dimensional annotations of the observation and variables are stored
    in the :attr:`obs` and :attr:`var` attributes as :class:`~pandas.DataFrame`\\ s.
    This is intended for metrics calculated over their axes.
    Multi-dimensional annotations are stored in :attr:`obsm` and :attr:`varm`,
    which are aligned to the objects observation and variable dimensions respectively.
    Square matrices representing graphs are stored in :attr:`obsp` and :attr:`varp`,
    with both of their own dimensions aligned to their associated axis.
    Additional measurements across both observations and variables are stored in
    :attr:`layers`.

    Indexing into an AnnData object can be performed by relative position
    with numeric indices (like pandas’ :meth:`~pandas.DataFrame.iloc`),
    or by labels (like :meth:`~pandas.DataFrame.loc`).
    To avoid ambiguity with numeric indexing into observations or variables,
    indexes of the AnnData object are converted to strings by the constructor.

    Subsetting an AnnData object by indexing into it will also subset its elements
    according to the dimensions they were aligned to.
    This means an operation like `adata[list_of_obs, :]` will also subset :attr:`obs`,
    :attr:`obsm`, and :attr:`layers`.

    Subsetting an AnnData object returns a view into the original object,
    meaning very little additional memory is used upon subsetting.
    This is achieved lazily, meaning that the constituent arrays are subset on access.
    Copying a view causes an equivalent “real” AnnData object to be generated.
    Attempting to modify a view (at any attribute except X) is handled
    in a copy-on-modify manner, meaning the object is initialized in place.
    Here’s an example::

        batch1 = adata[adata.obs["batch"] == "batch1", :]
        batch1.obs["value"] = 0  # This makes batch1 a “real” AnnData object

    At the end of this snippet: `adata` was not modified,
    and `batch1` is its own AnnData object with its own data.

    Similar to Bioconductor’s `ExpressionSet` and :mod:`scipy.sparse` matrices,
    subsetting an AnnData object retains the dimensionality of its constituent arrays.
    Therefore, unlike with the classes exposed by :mod:`pandas`, :mod:`numpy`,
    and `xarray`, there is no concept of a one dimensional AnnData object.
    AnnDatas always have two inherent dimensions, :attr:`obs` and :attr:`var`.
    Additionally, maintaining the dimensionality of the AnnData object allows for
    consistent handling of :mod:`scipy.sparse` matrices and :mod:`numpy` arrays.

    .. _statsmodels: http://www.statsmodels.org/stable/index.html
    .. _scikit-learn: http://scikit-learn.org/
    """

    _BACKED_ATTRS: ClassVar[list[str]] = ["X", "raw.X"]

    # backwards compat
    _H5_ALIASES: ClassVar[dict[str, set[str]]] = dict(
        X={"X", "_X", "data", "_data"},
        obs={"obs", "_obs", "smp", "_smp"},
        var={"var", "_var"},
        uns={"uns"},
        obsm={"obsm", "_obsm", "smpm", "_smpm"},
        varm={"varm", "_varm"},
        layers={"layers", "_layers"},
    )

    _H5_ALIASES_NAMES: ClassVar[dict[str, set[str]]] = dict(
        obs={"obs_names", "smp_names", "row_names", "index"},
        var={"var_names", "col_names", "index"},
    )

    _accessors: ClassVar[set[str]] = set()

    @old_positionals(
        "obsm",
        "varm",
        "layers",
        "raw",
        "dtype",
        "shape",
        "filename",
        "filemode",
        "asview",
    )
    def __init__(  # noqa: PLR0913
        self,
        X: XDataType | pd.DataFrame | None = None,
        obs: pd.DataFrame | Mapping[str, Iterable[Any]] | None = None,
        var: pd.DataFrame | Mapping[str, Iterable[Any]] | None = None,
        uns: Mapping[str, Any] | None = None,
        *,
        obsm: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
        varm: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
        layers: Mapping[str, XDataType] | None = None,
        raw: Mapping[str, Any] | None = None,
        dtype: np.dtype | type | str | None = None,
        shape: tuple[int, int] | None = None,
        filename: PathLike[str] | str | None = None,
        filemode: Literal["r", "r+"] | None = None,
        asview: bool = False,
        obsp: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
        varp: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
        oidx: Index1D | None = None,
        vidx: Index1D | None = None,
    ):
        # check for any multi-indices that aren’t later checked in coerce_array
        for attr, key in [(obs, "obs"), (var, "var"), (X, "X")]:
            if isinstance(attr, pd.DataFrame):
                raise_value_error_if_multiindex_columns(attr, key)
        if asview:
            if not isinstance(X, AnnData):
                msg = "`X` has to be an AnnData object."
                raise ValueError(msg)
            self._init_as_view(X, oidx, vidx)
        else:
            self._init_as_actual(
                X=X,
                obs=obs,
                var=var,
                uns=uns,
                obsm=obsm,
                varm=varm,
                raw=raw,
                layers=layers,
                dtype=dtype,
                shape=shape,
                obsp=obsp,
                varp=varp,
                filename=filename,
                filemode=filemode,
            )

    def _init_as_view(self, adata_ref: AnnData, oidx: Index, vidx: Index):
        if adata_ref.isbacked and adata_ref.is_view:
            msg = (
                "Currently, you cannot index repeatedly into a backed AnnData, "
                "that is, you cannot make a view of a view."
            )
            raise ValueError(msg)
        self._is_view = True
        if isinstance(oidx, int | np.integer):
            if not (-adata_ref.n_obs <= oidx < adata_ref.n_obs):
                msg = f"Observation index `{oidx}` is out of range."
                raise IndexError(msg)
            oidx += adata_ref.n_obs * (oidx < 0)
            oidx = slice(oidx, oidx + 1, 1)
        if isinstance(vidx, int | np.integer):
            if not (-adata_ref.n_vars <= vidx < adata_ref.n_vars):
                msg = f"Variable index `{vidx}` is out of range."
                raise IndexError(msg)
            vidx += adata_ref.n_vars * (vidx < 0)
            vidx = slice(vidx, vidx + 1, 1)
        if adata_ref.is_view:
            prev_oidx, prev_vidx = adata_ref._oidx, adata_ref._vidx
            adata_ref = adata_ref._adata_ref
            oidx, vidx = _resolve_idxs((prev_oidx, prev_vidx), (oidx, vidx), adata_ref)
        # self._adata_ref is never a view
        self._adata_ref = adata_ref
        self._oidx = oidx
        self._vidx = vidx
        # the file is the same as of the reference object
        self.file = adata_ref.file
        # views on attributes of adata_ref
        obs_sub = adata_ref.obs.iloc[oidx]
        var_sub = adata_ref.var.iloc[vidx]
        # fix categories
        uns = copy(adata_ref._uns)
        if settings.remove_unused_categories:
            self._remove_unused_categories(adata_ref.obs, obs_sub, uns)
            self._remove_unused_categories(adata_ref.var, var_sub, uns)
        # set attributes
        self._obs = as_view(obs_sub, view_args=(self, "obs"))
        self._var = as_view(var_sub, view_args=(self, "var"))
        self._uns = uns

        # set data
        if self.isbacked:
            self._X = None

        # set raw, easy, as it’s immutable anyways...
        if adata_ref._raw is not None:
            # slicing along variables axis is ignored
            self._raw = adata_ref.raw[oidx]
            self._raw._adata = self
        else:
            self._raw = None

    def _init_as_actual(  # noqa: PLR0912, PLR0913, PLR0915
        self,
        X=None,
        *,
        obs=None,
        var=None,
        uns=None,
        obsm=None,
        varm=None,
        varp=None,
        obsp=None,
        raw=None,
        layers=None,
        dtype=None,
        shape=None,
        filename=None,
        filemode=None,
    ):
        # view attributes
        self._is_view = False
        self._adata_ref = None
        self._oidx = None
        self._vidx = None

        # ----------------------------------------------------------------------
        # various ways of initializing the data
        # ----------------------------------------------------------------------

        # If X is a data frame, we store its indices for verification
        x_indices = []

        # init from file
        if filename is not None:
            self.file = AnnDataFileManager(self, filename, filemode)
        else:
            self.file = AnnDataFileManager(self, None)

            # init from AnnData
            if isinstance(X, AnnData):
                if any((obs, var, uns, obsm, varm, obsp, varp)):
                    msg = "If `X` is a dict no further arguments must be provided."
                    raise ValueError(msg)
                X, obs, var, uns, obsm, varm, obsp, varp, layers, raw = (
                    X._X,
                    X.obs,
                    X.var,
                    X.uns,
                    X.obsm,
                    X.varm,
                    X.obsp,
                    X.varp,
                    X.layers,
                    X.raw,
                )

            # init from DataFrame
            elif isinstance(X, pd.DataFrame):
                # to verify index matching, we wait until obs and var are DataFrames
                if obs is None:
                    obs = pd.DataFrame(index=X.index)
                elif not isinstance(X.index, pd.RangeIndex):
                    x_indices.append(("obs", "index", X.index.astype(str)))
                if var is None:
                    var = pd.DataFrame(index=X.columns)
                elif not isinstance(X.columns, pd.RangeIndex):
                    x_indices.append(("var", "columns", X.columns.astype(str)))
                X = ensure_df_homogeneous(X, "X")

        # ----------------------------------------------------------------------
        # actually process the data
        # ----------------------------------------------------------------------

        # check data type of X
        if X is not None:
            X = coerce_array(X, name="X")
            if shape is not None:
                msg = "`shape` needs to be `None` if `X` is not `None`."
                raise ValueError(msg)
            _check_2d_shape(X)
            # if type doesn’t match, a copy is made, otherwise, use a view
            if dtype is not None:
                msg = (
                    "The dtype argument is deprecated and will be removed in late 2024."
                )
                warnings.warn(msg, FutureWarning, stacklevel=3)
                if issparse(X) or isinstance(X, ma.MaskedArray):
                    # TODO: maybe use view on data attribute of sparse matrix
                    #       as in readwrite.read_10x_h5
                    if X.dtype != np.dtype(dtype):
                        X = X.astype(dtype)
                elif isinstance(X, ZarrArray | DaskArray):
                    X = X.astype(dtype)
                else:  # is np.ndarray or a subclass, convert to true np.ndarray
                    X = np.asarray(X, dtype)
            # data matrix and shape
            self._X = X
            n_obs, n_vars = X.shape
            source = "X"
        else:
            self._X = None
            n_obs, n_vars = (
                shape
                if shape is not None
                else _infer_shape(
                    obs, var, obsm=obsm, varm=varm, layers=layers, obsp=obsp, varp=varp
                )
            )
            source = "shape"

        # annotations
        self._obs = _gen_dataframe(
            obs, ["obs_names", "row_names"], source=source, attr="obs", length=n_obs
        )
        self._var = _gen_dataframe(
            var, ["var_names", "col_names"], source=source, attr="var", length=n_vars
        )

        # now we can verify if indices match!
        for attr_name, x_name, idx in x_indices:
            attr = getattr(self, attr_name)
            if isinstance(attr.index, pd.RangeIndex):
                attr.index = idx
            elif not idx.equals(attr.index):
                msg = f"Index of {attr_name} must match {x_name} of X."
                raise ValueError(msg)

        # unstructured annotations
        self.uns = uns or OrderedDict()

        self.obsm = obsm
        self.varm = varm

        self.obsp = obsp
        self.varp = varp

        # Backwards compat for connectivities matrices in uns["neighbors"]
        _move_adj_mtx({"uns": self._uns, "obsp": self._obsp})
        self._check_dimensions()
        if settings.check_uniqueness:
            self._check_uniqueness()

        if self.filename:
            assert not isinstance(raw, Raw), (
                "got raw from other adata but also filename?"
            )
            if {"raw", "raw.X"} & set(self.file):
                raw = dict(X=None, **raw)
        if not raw:
            self._raw = None
        elif isinstance(raw, Mapping):
            self._raw = Raw(self, **raw)
        else:  # is a Raw from another AnnData
            self._raw = Raw(self, raw._X, raw.var, raw.varm)

        # clean up old formats
        self._clean_up_old_format(uns)

        # layers
        self.layers = layers

    @old_positionals("show_stratified", "with_disk")
    def __sizeof__(
        self, *, show_stratified: bool = False, with_disk: bool = False
    ) -> int:
        def get_size(X) -> int:
            def cs_to_bytes(X) -> int:
                return int(X.data.nbytes + X.indptr.nbytes + X.indices.nbytes)

            if isinstance(X, h5py.Dataset) and with_disk:
                return int(np.array(X.shape).prod() * X.dtype.itemsize)
            elif isinstance(X, BaseCompressedSparseDataset) and with_disk:
                return cs_to_bytes(X._to_backed())
            elif issparse(X):
                return cs_to_bytes(X)
            else:
                return X.__sizeof__()

        sizes = {}
        attrs = ["X", "_obs", "_var"]
        attrs_multi = ["_uns", "_obsm", "_varm", "varp", "_obsp", "_layers"]
        for attr in attrs + attrs_multi:
            if attr in attrs_multi:
                keys = getattr(self, attr).keys()
                s = sum(get_size(getattr(self, attr)[k]) for k in keys)
            else:
                s = get_size(getattr(self, attr))
            if s > 0 and show_stratified:
                from tqdm import tqdm

                print(
                    f"Size of {attr.replace('_', '.'):<7}: {tqdm.format_sizeof(s, 'B')}"
                )
            sizes[attr] = s
        return sum(sizes.values())

    def _gen_repr(self, n_obs, n_vars) -> str:
        backed_at = f" backed at {str(self.filename)!r}" if self.isbacked else ""
        descr = f"AnnData object with n_obs × n_vars = {n_obs} × {n_vars}{backed_at}"
        for attr in [
            "obs",
            "var",
            "uns",
            "obsm",
            "varm",
            "layers",
            "obsp",
            "varp",
        ]:
            keys = getattr(self, attr).keys()
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return descr

    def __repr__(self) -> str:
        if self.is_view:
            return "View of " + self._gen_repr(self.n_obs, self.n_vars)
        else:
            return self._gen_repr(self.n_obs, self.n_vars)

    def __eq__(self, other):
        """Equality testing"""
        msg = (
            "Equality comparisons are not supported for AnnData objects, "
            "instead compare the desired attributes."
        )
        raise NotImplementedError(msg)

    @property
    def shape(self) -> tuple[int, int]:
        """Shape of data matrix (:attr:`n_obs`, :attr:`n_vars`)."""
        return self.n_obs, self.n_vars

    @property
    def X(self) -> XDataType | None:
        """Data matrix of shape :attr:`n_obs` × :attr:`n_vars`."""
        if self.isbacked:
            if not self.file.is_open:
                self.file.open()
            X = self.file["X"]
            if isinstance(X, h5py.Group):
                X = sparse_dataset(X)
            # This is so that we can index into a backed dense dataset with
            # indices that aren’t strictly increasing
            if self.is_view:
                X = _subset(X, (self._oidx, self._vidx))
        elif self.is_view and self._adata_ref.X is None:
            X = None
        elif self.is_view:
            X = as_view(
                _subset(self._adata_ref.X, (self._oidx, self._vidx)),
                ElementRef(self, "X"),
            )
        else:
            X = self._X
        return X
        # if self.n_obs == 1 and self.n_vars == 1:
        #     return X[0, 0]
        # elif self.n_obs == 1 or self.n_vars == 1:
        #     if issparse(X): X = X.toarray()
        #     return X.flatten()
        # else:
        #     return X

    @X.setter
    def X(self, value: XDataType | None):  # noqa: PLR0912
        if value is None:
            if self.isbacked:
                msg = "Cannot currently remove data matrix from backed object."
                raise NotImplementedError(msg)
            if self.is_view:
                self._init_as_actual(self.copy())
            self._X = None
            return
        value = coerce_array(value, name="X", allow_array_like=True)

        # If indices are both arrays, we need to modify them
        # so we don’t set values like coordinates
        # This can occur if there are successive views
        if (
            self.is_view
            and isinstance(self._oidx, np.ndarray)
            and isinstance(self._vidx, np.ndarray)
        ):
            oidx, vidx = np.ix_(self._oidx, self._vidx)
        else:
            oidx, vidx = self._oidx, self._vidx
        if (
            np.isscalar(value)
            or (hasattr(value, "shape") and (self.shape == value.shape))
            or (self.n_vars == 1 and self.n_obs == len(value))
            or (self.n_obs == 1 and self.n_vars == len(value))
        ):
            if not np.isscalar(value):
                if self.is_view and any(
                    isinstance(idx, np.ndarray)
                    and len(np.unique(idx)) != len(idx.ravel())
                    for idx in [oidx, vidx]
                ):
                    msg = (
                        "You are attempting to set `X` to a matrix on a view which has non-unique indices. "
                        "The resulting `adata.X` will likely not equal the value to which you set it. "
                        "To avoid this potential issue, please make a copy of the data first. "
                        "In the future, this operation will throw an error."
                    )
                    warnings.warn(msg, FutureWarning, stacklevel=1)
                if self.shape != value.shape:
                    # For assigning vector of values to 2d array or matrix
                    # Not necessary for row of 2d array
                    value = value.reshape(self.shape)
            if self.isbacked:
                if self.is_view:
                    X = self.file["X"]
                    if isinstance(X, h5py.Group):
                        X = sparse_dataset(X)
                    X[oidx, vidx] = value
                else:
                    self._set_backed("X", value)
            elif self.is_view:
                if sparse.issparse(self._adata_ref._X) and isinstance(
                    value, np.ndarray
                ):
                    if isinstance(self._adata_ref.X, CSArray):
                        memory_class = sparse.coo_array
                    else:
                        memory_class = sparse.coo_matrix
                    value = memory_class(value)
                elif sparse.issparse(value) and isinstance(
                    self._adata_ref._X, np.ndarray
                ):
                    warnings.warn(
                        "Trying to set a dense array with a sparse array on a view."
                        "Densifying the sparse array."
                        "This may incur excessive memory usage",
                        stacklevel=2,
                    )
                    value = value.toarray()
                warnings.warn(
                    "Modifying `X` on a view results in data being overridden",
                    ImplicitModificationWarning,
                    stacklevel=2,
                )
                self._adata_ref._X[oidx, vidx] = value
            else:
                self._X = value
        else:
            msg = f"Data matrix has wrong shape {value.shape}, need to be {self.shape}."
            raise ValueError(msg)

    @X.deleter
    def X(self):
        self.X = None

    layers: AlignedMappingProperty[Layers | LayersView] = AlignedMappingProperty(
        "layers", Layers
    )
    """\
    Dictionary-like object with values of the same dimensions as :attr:`X`.

    Layers in AnnData are inspired by loompy’s :ref:`loomlayers`.

    Return the layer named `"unspliced"`::

        adata.layers["unspliced"]

    Create or replace the `"spliced"` layer::

        adata.layers["spliced"] = ...

    Assign the 10th column of layer `"spliced"` to the variable a::

        a = adata.layers["spliced"][:, 10]

    Delete the `"spliced"` layer::

        del adata.layers["spliced"]

    Return layers’ names::

        adata.layers.keys()
    """

    @property
    def raw(self) -> Raw:
        """\
        Store raw version of :attr:`X` and :attr:`var` as `.raw.X` and `.raw.var`.

        The :attr:`raw` attribute is initialized with the current content
        of an object by setting::

            adata.raw = adata.copy()

        Its content can be deleted::

            adata.raw = None
            # or
            del adata.raw

        Upon slicing an AnnData object along the obs (row) axis, :attr:`raw`
        is also sliced. Slicing an AnnData object along the vars (columns) axis
        leaves :attr:`raw` unaffected. Note that you can call::

             adata.raw[:, 'orig_variable_name'].X

        to retrieve the data associated with a variable that might have been
        filtered out or "compressed away" in :attr:`X`.
        """
        return self._raw

    @raw.setter
    def raw(self, value: AnnData):
        if value is None:
            del self.raw
        elif not isinstance(value, AnnData):
            msg = "Can only init raw attribute with an AnnData object."
            raise ValueError(msg)
        else:
            if self.is_view:
                self._init_as_actual(self.copy())
            self._raw = Raw(self, X=value.X, var=value.var, varm=value.varm)

    @raw.deleter
    def raw(self):
        if self.is_view:
            self._init_as_actual(self.copy())
        self._raw = None

    @property
    def n_obs(self) -> int:
        """Number of observations."""
        return len(self.obs_names)

    @property
    def n_vars(self) -> int:
        """Number of variables/features."""
        return len(self.var_names)

    def _set_dim_df(self, value: pd.DataFrame | XDataset, attr: Literal["obs", "var"]):
        value = _gen_dataframe(
            value,
            [f"{attr}_names", f"{'row' if attr == 'obs' else 'col'}_names"],
            source="shape",
            attr=attr,
            length=self.n_obs if attr == "obs" else self.n_vars,
        )
        raise_value_error_if_multiindex_columns(value, attr)
        value_idx = self._prep_dim_index(value.index, attr)
        if self.is_view:
            self._init_as_actual(self.copy())
        setattr(self, f"_{attr}", value)
        self._set_dim_index(value_idx, attr)
        if not len(value.columns):
            value.columns = value.columns.astype(str)

    def _prep_dim_index(self, value, attr: str) -> pd.Index:
        """Prepares index to be uses as obs_names or var_names for AnnData object.AssertionError

        If a pd.Index is passed, this will use a reference, otherwise a new index object is created.
        """
        if self.shape[attr == "var"] != len(value):
            msg = f"Length of passed value for {attr}_names is {len(value)}, but this AnnData has shape: {self.shape}"
            raise ValueError(msg)
        if isinstance(value, pd.Index) and not isinstance(value.name, str | type(None)):
            msg = (
                f"AnnData expects .{attr}.index.name to be a string or None, "
                f"but you passed a name of type {type(value.name).__name__!r}"
            )
            raise ValueError(msg)
        else:
            value = pd.Index(value)
            if not isinstance(value.name, str | type(None)):
                value.name = None
        if (
            len(value) > 0
            and not isinstance(value, pd.RangeIndex)
            and infer_dtype(value) not in {"string", "bytes"}
        ):
            sample = list(value[: min(len(value), 5)])
            msg = dedent(
                f"""
                AnnData expects .{attr}.index to contain strings, but got values like:
                    {sample}

                    Inferred to be: {infer_dtype(value)}
                """
            )
            warnings.warn(msg, stacklevel=2)
        return value

    def _set_dim_index(self, value: pd.Index, attr: str):
        # Assumes _prep_dim_index has been run
        if self.is_view:
            self._init_as_actual(self.copy())
        getattr(self, attr).index = value
        for v in getattr(self, f"_{attr}m").values():
            if isinstance(v, pd.DataFrame):
                v.index = value

    @property
    def obs(self) -> pd.DataFrame | Dataset2D:
        """One-dimensional annotation of observations (`pd.DataFrame`)."""
        return self._obs

    @obs.setter
    def obs(self, value: pd.DataFrame | XDataset):
        self._set_dim_df(value, "obs")

    @obs.deleter
    def obs(self):
        self.obs = pd.DataFrame({}, index=self.obs_names)

    @property
    def obs_names(self) -> pd.Index:
        """Names of observations (alias for `.obs.index`)."""
        return self.obs.index

    @obs_names.setter
    def obs_names(self, names: Sequence[str]):
        names = self._prep_dim_index(names, "obs")
        self._set_dim_index(names, "obs")

    @property
    def var(self) -> pd.DataFrame | Dataset2D:
        """One-dimensional annotation of variables/ features (`pd.DataFrame`)."""
        return self._var

    @var.setter
    def var(self, value: pd.DataFrame | XDataset):
        self._set_dim_df(value, "var")

    @var.deleter
    def var(self):
        self.var = pd.DataFrame({}, index=self.var_names)

    @property
    def var_names(self) -> pd.Index:
        """Names of variables (alias for `.var.index`)."""
        return self.var.index

    @var_names.setter
    def var_names(self, names: Sequence[str]):
        names = self._prep_dim_index(names, "var")
        self._set_dim_index(names, "var")

    @property
    def uns(self) -> MutableMapping:
        """Unstructured annotation (ordered dictionary)."""
        uns = self._uns
        if self.is_view:
            uns = DictView(uns, view_args=(self, "_uns"))
        return uns

    @uns.setter
    def uns(self, value: MutableMapping):
        if not isinstance(value, MutableMapping):
            msg = "Only mutable mapping types (e.g. dict) are allowed for `.uns`."
            raise ValueError(msg)
        if isinstance(value, DictView):
            value = value.copy()
        if self.is_view:
            self._init_as_actual(self.copy())
        self._uns = value

    @uns.deleter
    def uns(self):
        self.uns = OrderedDict()

    obsm: AlignedMappingProperty[AxisArrays | AxisArraysView] = AlignedMappingProperty(
        "obsm", AxisArrays, 0
    )
    """\
    Multi-dimensional annotation of observations
    (mutable structured :class:`~numpy.ndarray`).

    Stores for each key a two or higher-dimensional :class:`~numpy.ndarray`
    of length `n_obs`.
    Is sliced with `data` and `obs` but behaves otherwise like a :term:`mapping`.
    """

    varm: AlignedMappingProperty[AxisArrays | AxisArraysView] = AlignedMappingProperty(
        "varm", AxisArrays, 1
    )
    """\
    Multi-dimensional annotation of variables/features
    (mutable structured :class:`~numpy.ndarray`).

    Stores for each key a two or higher-dimensional :class:`~numpy.ndarray`
    of length `n_vars`.
    Is sliced with `data` and `var` but behaves otherwise like a :term:`mapping`.
    """

    obsp: AlignedMappingProperty[PairwiseArrays | PairwiseArraysView] = (
        AlignedMappingProperty("obsp", PairwiseArrays, 0)
    )
    """\
    Pairwise annotation of observations,
    a mutable mapping with array-like values.

    Stores for each key a two or higher-dimensional :class:`~numpy.ndarray`
    whose first two dimensions are of length `n_obs`.
    Is sliced with `data` and `obs` but behaves otherwise like a :term:`mapping`.
    """

    varp: AlignedMappingProperty[PairwiseArrays | PairwiseArraysView] = (
        AlignedMappingProperty("varp", PairwiseArrays, 1)
    )
    """\
    Pairwise annotation of variables/features,
    a mutable mapping with array-like values.

    Stores for each key a two or higher-dimensional :class:`~numpy.ndarray`
    whose first two dimensions are of length `n_var`.
    Is sliced with `data` and `var` but behaves otherwise like a :term:`mapping`.
    """

    def obs_keys(self) -> list[str]:
        """List keys of observation annotation :attr:`obs`."""
        return self._obs.keys().tolist()

    def var_keys(self) -> list[str]:
        """List keys of variable annotation :attr:`var`."""
        return self._var.keys().tolist()

    def obsm_keys(self) -> list[str]:
        """List keys of observation annotation :attr:`obsm`."""
        return list(self.obsm.keys())

    def varm_keys(self) -> list[str]:
        """List keys of variable annotation :attr:`varm`."""
        return list(self.varm.keys())

    def uns_keys(self) -> list[str]:
        """List keys of unstructured annotation."""
        return sorted(self._uns.keys())

    @property
    def isbacked(self) -> bool:
        """`True` if object is backed on disk, `False` otherwise."""
        return self.filename is not None

    @property
    def is_view(self) -> bool:
        """`True` if object is view of another AnnData object, `False` otherwise."""
        return self._is_view

    @property
    def filename(self) -> Path | None:
        """\
        Change to backing mode by setting the filename of a `.h5ad` file.

        - Setting the filename writes the stored data to disk.
        - Setting the filename when the filename was previously another name
          moves the backing file from the previous file to the new file.
          If you want to copy the previous file, use `copy(filename='new_filename')`.
        """
        return self.file.filename

    @filename.setter
    def filename(self, filename: PathLike[str] | str | None):
        # convert early for later comparison
        filename = None if filename is None else Path(filename)
        # change from backing-mode back to full loading into memory
        if filename is None:
            if self.filename is not None:
                self.file._to_memory_mode()
            else:
                # both filename and self.filename are None
                # do nothing
                return
        else:
            if self.filename is not None:
                if self.filename != filename:
                    # write the content of self to the old file
                    # and close the file
                    self.write()
                    self.filename.rename(filename)
                else:
                    # do nothing
                    return
            else:
                # change from memory to backing-mode
                # write the content of self to disk
                as_dense = ("X", "raw/X") if self.raw is not None else ("X",)
                self.write(filename, as_dense=as_dense)
            # open new file for accessing
            self.file.open(filename, "r+")
            # as the data is stored on disk, we can safely set self._X to None
            self._X = None

    def _set_backed(self, attr, value):
        from .._io.utils import write_attribute

        write_attribute(self.file._file, attr, value)

    def _normalize_indices(self, index: Index | None) -> tuple[slice, slice]:
        return _normalize_indices(index, self.obs_names, self.var_names)

    # TODO: this is not quite complete...
    def __delitem__(self, index: Index):
        obs, var = self._normalize_indices(index)
        # TODO: does this really work?
        if not self.isbacked:
            del self._X[obs, var]
        else:
            X = self.file["X"]
            del X[obs, var]
            self._set_backed("X", X)
        if var == slice(None):
            del self._obs.iloc[obs, :]
        if obs == slice(None):
            del self._var.iloc[var, :]

    def __getitem__(self, index: Index) -> AnnData:
        """Returns a sliced view of the object."""
        oidx, vidx = self._normalize_indices(index)
        return AnnData(self, oidx=oidx, vidx=vidx, asview=True)

    @singledispatchmethod
    @staticmethod
    def _remove_unused_categories(
        df_full: pd.DataFrame, df_sub: pd.DataFrame, uns: dict[str, Any]
    ):
        for k in df_full:
            if not isinstance(df_full[k].dtype, pd.CategoricalDtype):
                continue
            all_categories = df_full[k].cat.categories
            with pd.option_context("mode.chained_assignment", None):
                df_sub[k] = df_sub[k].cat.remove_unused_categories()
            # also correct the colors...
            color_key = f"{k}_colors"
            if color_key not in uns:
                continue
            color_vec = uns[color_key]
            if np.array(color_vec).ndim == 0:
                # Make 0D arrays into 1D ones
                uns[color_key] = np.array(color_vec)[(None,)]
            elif len(color_vec) != len(all_categories):
                # Reset colors
                del uns[color_key]
            else:
                idx = np.where(np.isin(all_categories, df_sub[k].cat.categories))[0]
                uns[color_key] = np.array(color_vec)[(idx,)]

    def rename_categories(self, key: str, categories: Sequence[Any]):
        """\
        Rename categories of annotation `key` in :attr:`obs`, :attr:`var`,
        and :attr:`uns`.

        Only supports passing a list/array-like `categories` argument.

        Besides calling `self.obs[key].cat.categories = categories` –
        similar for :attr:`var` - this also renames categories in unstructured
        annotation that uses the categorical annotation `key`.

        Parameters
        ----------
        key
             Key for observations or variables annotation.
        categories
             New categories, the same number as the old categories.
        """
        if isinstance(categories, Mapping):
            msg = "Only list-like `categories` is supported."
            raise ValueError(msg)
        if key in self.obs:
            old_categories = self.obs[key].cat.categories.tolist()
            self.obs[key] = self.obs[key].cat.rename_categories(categories)
        elif key in self.var:
            old_categories = self.var[key].cat.categories.tolist()
            self.var[key] = self.var[key].cat.rename_categories(categories)
        else:
            msg = f"{key} is neither in `.obs` nor in `.var`."
            raise ValueError(msg)
        # this is not a good solution
        # but depends on the scanpy conventions for storing the categorical key
        # as `groupby` in the `params` slot
        for k1, v1 in self.uns.items():
            if not (
                isinstance(v1, Mapping)
                and "params" in v1
                and "groupby" in v1["params"]
                and v1["params"]["groupby"] == key
            ):
                continue
            for k2, v2 in v1.items():
                # picks out the recarrays that are named according to the old
                # categories
                if isinstance(v2, np.ndarray) and v2.dtype.names is not None:
                    if list(v2.dtype.names) == old_categories:
                        self.uns[k1][k2].dtype.names = categories
                    else:
                        logger.warning(
                            f"Omitting {k1}/{k2} as old categories do not match."
                        )

    def strings_to_categoricals(self, df: pd.DataFrame | None = None):
        """\
        Transform string annotations to categoricals.

        Only affects string annotations that lead to less categories than the
        total number of observations.

        Params
        ------
        df
            If `df` is `None`, modifies both :attr:`obs` and :attr:`var`,
            otherwise modifies `df` inplace.

        Notes
        -----
        Turns the view of an :class:`~anndata.AnnData` into an actual
        :class:`~anndata.AnnData`.
        """
        dont_modify = False  # only necessary for backed views
        if df is None:
            dfs = [self.obs, self.var]
            if self.is_view and self.isbacked:
                dont_modify = True
        else:
            dfs = [df]
        del df

        for df in dfs:
            string_cols = [
                key for key in df.columns if infer_dtype(df[key]) == "string"
            ]
            for key in string_cols:
                c = pd.Categorical(df[key])
                # TODO: We should only check if non-null values are unique, but
                # this would break cases where string columns with nulls could
                # be written as categorical, but not as string.
                # Possible solution: https://github.com/scverse/anndata/issues/504
                if len(c.categories) >= len(c):
                    continue
                # Ideally this could be done inplace
                sorted_categories = natsorted(c.categories)
                if not np.array_equal(c.categories, sorted_categories):
                    c = c.reorder_categories(sorted_categories)
                if dont_modify:
                    msg = (
                        "Please call `.strings_to_categoricals()` on full "
                        "AnnData, not on this view. You might encounter this"
                        "error message while copying or writing to disk."
                    )
                    raise RuntimeError(msg)
                df[key] = c
                logger.info(f"... storing {key!r} as categorical")

    _sanitize = strings_to_categoricals  # backwards compat

    def _inplace_subset_var(self, index: Index1D):
        """\
        Inplace subsetting along variables dimension.

        Same as `adata = adata[:, index]`, but inplace.
        """
        adata_subset = self[:, index].copy()

        self._init_as_actual(adata_subset)

    def _inplace_subset_obs(self, index: Index1D):
        """\
        Inplace subsetting along variables dimension.

        Same as `adata = adata[index, :]`, but inplace.
        """
        adata_subset = self[index].copy()

        self._init_as_actual(adata_subset)

    # TODO: Update, possibly remove
    def __setitem__(self, index: Index, val: float | XDataType):
        if self.is_view:
            msg = "Object is view and cannot be accessed with `[]`."
            raise ValueError(msg)
        obs, var = self._normalize_indices(index)
        if not self.isbacked:
            self._X[obs, var] = val
        else:
            X = self.file["X"]
            X[obs, var] = val
            self._set_backed("X", X)

    def __len__(self) -> int:
        return self.shape[0]

    def transpose(self) -> AnnData:
        """\
        Transpose whole object.

        Data matrix is transposed, observations and variables are interchanged.
        Ignores `.raw`.
        """
        from anndata.compat import _safe_transpose

        X = self.X if not self.isbacked else self.file["X"]
        if self.is_view:
            msg = (
                "You’re trying to transpose a view of an `AnnData`, "
                "which is currently not implemented. Call `.copy()` before transposing."
            )
            raise ValueError(msg)

        return AnnData(
            X=_safe_transpose(X) if X is not None else None,
            layers={k: _safe_transpose(v) for k, v in self.layers.items()},
            obs=self.var,
            var=self.obs,
            uns=self._uns,
            obsm=self.varm,
            varm=self.obsm,
            obsp=self.varp,
            varp=self.obsp,
            filename=self.filename,
        )

    T = property(transpose)

    def to_df(self, layer: str | None = None) -> pd.DataFrame:
        """\
        Generate shallow :class:`~pandas.DataFrame`.

        The data matrix :attr:`X` is returned as
        :class:`~pandas.DataFrame`, where :attr:`obs_names` initializes the
        index, and :attr:`var_names` the columns.

        * No annotations are maintained in the returned object.
        * The data matrix is densified in case it is sparse.

        Params
        ------
        layer
            Key for `.layers`.

        Returns
        -------
        Pandas DataFrame of specified data matrix.
        """
        if layer is not None:
            X = self.layers[layer]
        elif not self._has_X():
            msg = "X is None, cannot convert to dataframe."
            raise ValueError(msg)
        else:
            X = self.X
        if issparse(X):
            X = X.toarray()
        return pd.DataFrame(X, index=self.obs_names, columns=self.var_names)

    def _get_X(self, *, use_raw: bool = False, layer: str | None = None):
        """\
        Convenience method for getting expression values
        with common arguments and error handling.
        """
        is_layer = layer is not None
        if use_raw and is_layer:
            msg = (
                "Cannot use expression from both layer and raw. You provided:"
                f"`use_raw={use_raw}` and `layer={layer}`"
            )
            raise ValueError(msg)
        if is_layer:
            return self.layers[layer]
        elif use_raw:
            if self.raw is None:
                msg = "This AnnData doesn’t have a value in `.raw`."
                raise ValueError(msg)
            return self.raw.X
        else:
            return self.X

    def obs_vector(self, k: str, *, layer: str | None = None) -> np.ndarray:
        """\
        Convenience function for returning a 1 dimensional ndarray of values
        from :attr:`X`, :attr:`layers`\\ `[k]`, or :attr:`obs`.

        Made for convenience, not performance.
        Intentionally permissive about arguments, for easy iterative use.

        Params
        ------
        k
            Key to use. Should be in :attr:`var_names` or :attr:`obs`\\ `.columns`.
        layer
            What layer values should be returned from. If `None`, :attr:`X` is used.

        Returns
        -------
        A one dimensional ndarray, with values for each obs in the same order
        as :attr:`obs_names`.
        """
        if layer == "X":
            if "X" in self.layers:
                pass
            else:
                msg = (
                    "In a future version of AnnData, access to `.X` by passing"
                    " `layer='X'` will be removed. Instead pass `layer=None`."
                )
                warnings.warn(msg, FutureWarning, stacklevel=2)
                layer = None
        return get_vector(self, k, "obs", "var", layer=layer)

    def var_vector(self, k, *, layer: str | None = None) -> np.ndarray:
        """\
        Convenience function for returning a 1 dimensional ndarray of values
        from :attr:`X`, :attr:`layers`\\ `[k]`, or :attr:`obs`.

        Made for convenience, not performance. Intentionally permissive about
        arguments, for easy iterative use.

        Params
        ------
        k
            Key to use. Should be in :attr:`obs_names` or :attr:`var`\\ `.columns`.
        layer
            What layer values should be returned from. If `None`, :attr:`X` is used.

        Returns
        -------
        A one dimensional ndarray, with values for each var in the same order
        as :attr:`var_names`.
        """
        if layer == "X":
            if "X" in self.layers:
                pass
            else:
                msg = (
                    "In a future version of AnnData, access to `.X` by passing "
                    "`layer='X'` will be removed. Instead pass `layer=None`."
                )
                warnings.warn(msg, FutureWarning, stacklevel=2)
                layer = None
        return get_vector(self, k, "var", "obs", layer=layer)

    @deprecated("obs_vector")
    def _get_obs_array(self, k, use_raw=False, layer=None):  # noqa: FBT002
        """\
        Get an array from the layer (default layer='X') along the :attr:`obs`
        dimension by first looking up `obs.keys` and then :attr:`obs_names`.
        """
        if not use_raw or k in self.obs.columns:
            return self.obs_vector(k=k, layer=layer)
        else:
            return self.raw.obs_vector(k)

    @deprecated("var_vector")
    def _get_var_array(self, k, use_raw=False, layer=None):  # noqa: FBT002
        """\
        Get an array from the layer (default layer='X') along the :attr:`var`
        dimension by first looking up `var.keys` and then :attr:`var_names`.
        """
        if not use_raw or k in self.var.columns:
            return self.var_vector(k=k, layer=layer)
        else:
            return self.raw.var_vector(k)

    def _mutated_copy(self, **kwargs):
        """Creating AnnData with attributes optionally specified via kwargs."""
        if self.isbacked and (
            "X" not in kwargs or (self.raw is not None and "raw" not in kwargs)
        ):
            msg = (
                "This function does not currently handle backed objects "
                "internally, this should be dealt with before."
            )
            raise NotImplementedError(msg)
        new = {}

        for key in ["obs", "var", "obsm", "varm", "obsp", "varp", "layers"]:
            if key in kwargs:
                new[key] = kwargs[key]
            else:
                new[key] = getattr(self, key).copy()
        if "X" in kwargs:
            new["X"] = kwargs["X"]
        elif self._has_X():
            new["X"] = self.X.copy()
        if "uns" in kwargs:
            new["uns"] = kwargs["uns"]
        else:
            new["uns"] = deepcopy(self._uns)
        if "raw" in kwargs:
            new["raw"] = kwargs["raw"]
        elif self.raw is not None:
            new["raw"] = self.raw.copy()
        return AnnData(**new)

    @old_positionals("copy")
    def to_memory(self, *, copy: bool = False) -> AnnData:
        """Return a new AnnData object with all backed arrays loaded into memory.

        Params
        ------
        copy
            Whether the arrays that are already in-memory should be copied.

        Example
        -------

        .. code:: python

            import anndata
            backed = anndata.io.read_h5ad("file.h5ad", backed="r")
            mem = backed[backed.obs["cluster"] == "a", :].to_memory()
        """
        new = {}
        for attr_name in [
            "X",
            "obs",
            "var",
            "obsm",
            "varm",
            "obsp",
            "varp",
            "layers",
            "uns",
        ]:
            attr = getattr(self, attr_name, None)
            if attr is not None:
                new[attr_name] = to_memory(attr, copy=copy)

        if self.raw is not None:
            new["raw"] = {
                "X": to_memory(self.raw.X, copy=copy),
                "var": to_memory(self.raw.var, copy=copy),
                "varm": to_memory(self.raw.varm, copy=copy),
            }

        if self.isbacked:
            self.file.close()

        return AnnData(**new)

    def copy(self, filename: PathLike[str] | str | None = None) -> AnnData:
        """Full copy, optionally on disk."""
        if not self.isbacked:
            if self.is_view and self._has_X():
                # TODO: How do I unambiguously check if this is a copy?
                # Subsetting this way means we don’t have to have a view type
                # defined for the matrix, which is needed for some of the
                # current distributed backend. Specifically Dask.
                return self._mutated_copy(
                    X=_subset(self._adata_ref.X, (self._oidx, self._vidx)).copy()
                )
            else:
                return self._mutated_copy()
        else:
            from ..io import read_h5ad, write_h5ad

            if filename is None:
                msg = (
                    "To copy an AnnData object in backed mode, "
                    "pass a filename: `.copy(filename='myfilename.h5ad')`. "
                    "To load the object into memory, use `.to_memory()`."
                )
                raise ValueError(msg)
            mode = self.file._filemode
            write_h5ad(filename, self)
            return read_h5ad(filename, backed=mode)

    @deprecated(
        "anndata.concat",
        add_msg="See the tutorial for concat at: "
        "https://anndata.readthedocs.io/en/latest/concatenation.html",
        hide=False,
    )
    def concatenate(
        self,
        *adatas: AnnData,
        join: str = "inner",
        batch_key: str = "batch",
        batch_categories: Sequence[Any] | None = None,
        uns_merge: str | None = None,
        index_unique: str | None = "-",
        fill_value=None,
    ) -> AnnData:
        """\
        Concatenate along the observations axis.

        The :attr:`uns`, :attr:`varm` and :attr:`obsm` attributes are ignored.

        Currently, this works only in `'memory'` mode.

        .. note::

            For more flexible and efficient concatenation, see: :func:`~anndata.concat`.

        Parameters
        ----------
        adatas
            AnnData matrices to concatenate with. Each matrix is referred to as
            a “batch”.
        join
            Use intersection (`'inner'`) or union (`'outer'`) of variables.
        batch_key
            Add the batch annotation to :attr:`obs` using this key.
        batch_categories
            Use these as categories for the batch annotation. By default, use increasing numbers.
        uns_merge
            Strategy to use for merging entries of uns. These strategies are applied recusivley.
            Currently implemented strategies include:

            * `None`: The default. The concatenated object will just have an empty dict for `uns`.
            * `"same"`: Only entries which have the same value in all AnnData objects are kept.
            * `"unique"`: Only entries which have one unique value in all AnnData objects are kept.
            * `"first"`: The first non-missing value is used.
            * `"only"`: A value is included if only one of the AnnData objects has a value at this
              path.
        index_unique
            Make the index unique by joining the existing index names with the
            batch category, using `index_unique='-'`, for instance. Provide
            `None` to keep existing indices.
        fill_value
            Scalar value to fill newly missing values in arrays with. Note: only applies to arrays
            and sparse matrices (not dataframes) and will only be used if `join="outer"`.

            .. note::
                If not provided, the default value is `0` for sparse matrices and `np.nan`
                for numpy arrays. See the examples below for more information.

        Returns
        -------
        :class:`~anndata.AnnData`
            The concatenated :class:`~anndata.AnnData`, where `adata.obs[batch_key]`
            stores a categorical variable labeling the batch.

        Notes
        -----

        .. warning::

           If you use `join='outer'` this fills 0s for sparse data when
           variables are absent in a batch. Use this with care. Dense data is
           filled with `NaN`. See the examples.

        Examples
        --------
        Joining on intersection of variables.

        >>> adata1 = AnnData(
        ...     np.array([[1, 2, 3], [4, 5, 6]]),
        ...     dict(obs_names=['s1', 's2'], anno1=['c1', 'c2']),
        ...     dict(var_names=['a', 'b', 'c'], annoA=[0, 1, 2]),
        ... )
        >>> adata2 = AnnData(
        ...     np.array([[1, 2, 3], [4, 5, 6]]),
        ...     dict(obs_names=['s3', 's4'], anno1=['c3', 'c4']),
        ...     dict(var_names=['d', 'c', 'b'], annoA=[0, 1, 2]),
        ... )
        >>> adata3 = AnnData(
        ...     np.array([[1, 2, 3], [4, 5, 6]]),
        ...     dict(obs_names=['s1', 's2'], anno2=['d3', 'd4']),
        ...     dict(var_names=['d', 'c', 'b'], annoA=[0, 2, 3], annoB=[0, 1, 2]),
        ... )
        >>> adata = adata1.concatenate(adata2, adata3)
        >>> adata
        AnnData object with n_obs × n_vars = 6 × 2
            obs: 'anno1', 'anno2', 'batch'
            var: 'annoA-0', 'annoA-1', 'annoA-2', 'annoB-2'
        >>> adata.X
        array([[2, 3],
               [5, 6],
               [3, 2],
               [6, 5],
               [3, 2],
               [6, 5]])
        >>> adata.obs
             anno1 anno2 batch
        s1-0    c1   NaN     0
        s2-0    c2   NaN     0
        s3-1    c3   NaN     1
        s4-1    c4   NaN     1
        s1-2   NaN    d3     2
        s2-2   NaN    d4     2
        >>> adata.var.T
                 b  c
        annoA-0  1  2
        annoA-1  2  1
        annoA-2  3  2
        annoB-2  2  1

        Joining on the union of variables.

        >>> outer = adata1.concatenate(adata2, adata3, join='outer')
        >>> outer
        AnnData object with n_obs × n_vars = 6 × 4
            obs: 'anno1', 'anno2', 'batch'
            var: 'annoA-0', 'annoA-1', 'annoA-2', 'annoB-2'
        >>> outer.var.T
                   a    b    c    d
        annoA-0  0.0  1.0  2.0  NaN
        annoA-1  NaN  2.0  1.0  0.0
        annoA-2  NaN  3.0  2.0  0.0
        annoB-2  NaN  2.0  1.0  0.0
        >>> outer.var_names
        Index(['a', 'b', 'c', 'd'], dtype='object')
        >>> outer.X
        array([[ 1.,  2.,  3., nan],
               [ 4.,  5.,  6., nan],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.]])
        >>> outer.X.sum(axis=0)
        array([nan, 25., 23., nan])
        >>> import pandas as pd
        >>> Xdf = pd.DataFrame(outer.X, columns=outer.var_names)
        >>> Xdf
             a    b    c    d
        0  1.0  2.0  3.0  NaN
        1  4.0  5.0  6.0  NaN
        2  NaN  3.0  2.0  1.0
        3  NaN  6.0  5.0  4.0
        4  NaN  3.0  2.0  1.0
        5  NaN  6.0  5.0  4.0
        >>> Xdf.sum()
        a     5.0
        b    25.0
        c    23.0
        d    10.0
        dtype: float64

        One way to deal with missing values is to use masked arrays:

        >>> from numpy import ma
        >>> outer.X = ma.masked_invalid(outer.X)
        >>> outer.X
        masked_array(
          data=[[1.0, 2.0, 3.0, --],
                [4.0, 5.0, 6.0, --],
                [--, 3.0, 2.0, 1.0],
                [--, 6.0, 5.0, 4.0],
                [--, 3.0, 2.0, 1.0],
                [--, 6.0, 5.0, 4.0]],
          mask=[[False, False, False,  True],
                [False, False, False,  True],
                [ True, False, False, False],
                [ True, False, False, False],
                [ True, False, False, False],
                [ True, False, False, False]],
          fill_value=1e+20)
        >>> outer.X.sum(axis=0).data
        array([ 5., 25., 23., 10.])

        The masked array is not saved but has to be reinstantiated after saving.

        >>> outer.write('./test.h5ad')
        >>> from anndata import read_h5ad
        >>> outer = read_h5ad('./test.h5ad')
        >>> outer.X
        array([[ 1.,  2.,  3., nan],
               [ 4.,  5.,  6., nan],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.],
               [nan,  3.,  2.,  1.],
               [nan,  6.,  5.,  4.]])

        For sparse data, everything behaves similarly,
        except that for `join='outer'`, zeros are added.

        >>> from scipy.sparse import csr_matrix
        >>> adata1 = AnnData(
        ...     csr_matrix([[0, 2, 3], [0, 5, 6]], dtype=np.float32),
        ...     dict(obs_names=['s1', 's2'], anno1=['c1', 'c2']),
        ...     dict(var_names=['a', 'b', 'c']),
        ... )
        >>> adata2 = AnnData(
        ...     csr_matrix([[0, 2, 3], [0, 5, 6]], dtype=np.float32),
        ...     dict(obs_names=['s3', 's4'], anno1=['c3', 'c4']),
        ...     dict(var_names=['d', 'c', 'b']),
        ... )
        >>> adata3 = AnnData(
        ... csr_matrix([[1, 2, 0], [0, 5, 6]], dtype=np.float32),
        ...     dict(obs_names=['s5', 's6'], anno2=['d3', 'd4']),
        ...     dict(var_names=['d', 'c', 'b']),
        ... )
        >>> adata = adata1.concatenate(adata2, adata3, join='outer')
        >>> adata.var_names
        Index(['a', 'b', 'c', 'd'], dtype='object')
        >>> adata.X.toarray()
        array([[0., 2., 3., 0.],
               [0., 5., 6., 0.],
               [0., 3., 2., 0.],
               [0., 6., 5., 0.],
               [0., 0., 2., 1.],
               [0., 6., 5., 0.]], dtype=float32)
        """
        from .merge import concat, merge_dataframes, merge_outer, merge_same

        if self.isbacked:
            msg = "Currently, concatenate only works in memory mode."
            raise ValueError(msg)

        if len(adatas) == 0:
            return self.copy()
        elif len(adatas) == 1 and not isinstance(adatas[0], AnnData):
            adatas = adatas[0]  # backwards compatibility
        all_adatas = (self, *adatas)

        out = concat(
            all_adatas,
            axis=0,
            join=join,
            label=batch_key,
            keys=batch_categories,
            uns_merge=uns_merge,
            fill_value=fill_value,
            index_unique=index_unique,
            pairwise=False,
        )

        # Backwards compat (some of this could be more efficient)
        # obs used to always be an outer join
        sparse_class = sparse.csr_matrix
        if any(isinstance(a.X, CSArray) for a in all_adatas):
            sparse_class = sparse.csr_array
        out.obs = concat(
            [AnnData(sparse_class(a.shape), obs=a.obs) for a in all_adatas],
            axis=0,
            join="outer",
            label=batch_key,
            keys=batch_categories,
            index_unique=index_unique,
        ).obs
        # Removing varm
        del out.varm
        # Implementing old-style merging of var
        if batch_categories is None:
            batch_categories = np.arange(len(all_adatas)).astype(str)
        pat = rf"-({'|'.join(batch_categories)})$"
        out.var = merge_dataframes(
            [a.var for a in all_adatas],
            out.var_names,
            partial(merge_outer, batch_keys=batch_categories, merge=merge_same),
        )
        out.var = out.var.iloc[
            :,
            (
                out.var.columns.str.extract(pat, expand=False)
                .fillna("")
                .argsort(kind="stable")
            ),
        ]

        return out

    def var_names_make_unique(self, join: str = "-"):
        # Important to go through the setter so obsm dataframes are updated too
        self.var_names = utils.make_index_unique(self.var.index, join)

    var_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    def obs_names_make_unique(self, join: str = "-"):
        # Important to go through the setter so obsm dataframes are updated too
        self.obs_names = utils.make_index_unique(self.obs.index, join)

    obs_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    def _check_uniqueness(self):
        if not self.obs.index.is_unique:
            utils.warn_names_duplicates("obs")
        if not self.var.index.is_unique:
            utils.warn_names_duplicates("var")

    def __contains__(self, key: Any):
        msg = "AnnData has no attribute __contains__, don’t check `in adata`."
        raise AttributeError(msg)

    def _check_dimensions(self, key=None):
        key = {"obsm", "varm"} if key is None else {key}
        if "obsm" in key and (
            not all(axis_len(o, 0) == self.n_obs for o in self.obsm.values())
            and len(self.obsm.dim_names) != self.n_obs
        ):
            msg = (
                "Observations annot. `obsm` must have number of rows of `X`"
                f" ({self.n_obs}), but has {len(self.obsm)} rows."
            )
            raise ValueError(msg)
        if "varm" in key and (
            not all(axis_len(v, 0) == self.n_vars for v in self.varm.values())
            and len(self.varm.dim_names) != self.n_vars
        ):
            msg = (
                "Variables annot. `varm` must have number of columns of `X`"
                f" ({self.n_vars}), but has {len(self.varm)} rows."
            )
            raise ValueError(msg)

    @old_positionals("compression", "compression_opts", "as_dense")
    def write_h5ad(
        self,
        filename: PathLike[str] | str | None = None,
        *,
        convert_strings_to_categoricals: bool = True,
        compression: Literal["gzip", "lzf"] | None = None,
        compression_opts: int | Any = None,
        as_dense: Sequence[str] = (),
    ):
        """\
        Write `.h5ad`-formatted hdf5 file.

        .. note::
           Setting compression to `'gzip'` can save disk space
           but will slow down writing and subsequent reading.
           Prior to v0.6.16, this was the default for parameter `compression`.

        Generally, if you have sparse data that are stored as a dense matrix,
        you can dramatically improve performance and reduce disk space
        by converting to a :class:`~scipy.sparse.csr_matrix`::

            from scipy.sparse import csr_matrix
            adata.X = csr_matrix(adata.X)

        Parameters
        ----------
        filename
            Filename of data file. Defaults to backing file.
        convert_strings_to_categoricals
            Convert string columns to categorical.
        compression
            For [`lzf`, `gzip`], see the h5py :ref:`dataset_compression`.

            Alternative compression filters such as `zstd` can be passed
            from the :doc:`hdf5plugin <hdf5plugin:usage>` library.
            Experimental.

            Usage example::

                import hdf5plugin
                adata.write_h5ad(
                    filename,
                    compression=hdf5plugin.FILTERS["zstd"]
                )

            .. note::
                Datasets written with hdf5plugin-provided compressors
                cannot be opened without first loading the hdf5plugin
                library using `import hdf5plugin`. When using alternative
                compression filters such as `zstd`, consider writing to
                `zarr` format instead of `h5ad`, as the `zarr` library
                provides a more transparent compression pipeline.

        compression_opts
            For [`lzf`, `gzip`], see the h5py :ref:`dataset_compression`.

            Alternative compression filters such as `zstd` can be configured
            using helpers from the :doc:`hdf5plugin <hdf5plugin:usage>`
            library. Experimental.

            Usage example (setting `zstd` compression level to 5)::

                import hdf5plugin
                adata.write_h5ad(
                    filename,
                    compression=hdf5plugin.FILTERS["zstd"],
                    compression_opts=hdf5plugin.Zstd(clevel=5).filter_options
                )

        as_dense
            Sparse arrays in AnnData object to write as dense. Currently only
            supports `X` and `raw/X`.
        """
        from ..io import write_h5ad

        if filename is None and not self.isbacked:
            msg = "Provide a filename!"
            raise ValueError(msg)
        if filename is None:
            filename = self.filename

        write_h5ad(
            Path(filename),
            self,
            convert_strings_to_categoricals=convert_strings_to_categoricals,
            compression=compression,
            compression_opts=compression_opts,
            as_dense=as_dense,
        )

        if self.isbacked:
            self.file.filename = filename

    write = write_h5ad  # a shortcut and backwards compat

    @old_positionals("skip_data", "sep")
    def write_csvs(
        self, dirname: PathLike[str] | str, *, skip_data: bool = True, sep: str = ","
    ):
        """\
        Write annotation to `.csv` files.

        It is not possible to recover the full :class:`~anndata.AnnData` from
        these files. Use :meth:`write` for this.

        Parameters
        ----------
        dirname
            Name of directory to which to export.
        skip_data
             Skip the data matrix :attr:`X`.
        sep
             Separator for the data.
        """
        from ..io import write_csvs

        write_csvs(dirname, self, skip_data=skip_data, sep=sep)

    @old_positionals("write_obsm_varm")
    def write_loom(
        self, filename: PathLike[str] | str, *, write_obsm_varm: bool = False
    ):
        """\
        Write `.loom`-formatted hdf5 file.

        Parameters
        ----------
        filename
            The filename.
        """
        from ..io import write_loom

        write_loom(filename, self, write_obsm_varm=write_obsm_varm)

    @old_positionals("chunks")
    def write_zarr(
        self,
        store: StoreLike,
        *,
        chunks: tuple[int, ...] | None = None,
        convert_strings_to_categoricals: bool = True,
    ):
        """\
        Write a hierarchical Zarr array store.

        Parameters
        ----------
        store
            The filename, a :class:`~typing.MutableMapping`, or a Zarr storage class.
        chunks
            Chunk shape.
        convert_strings_to_categoricals
            Convert string columns to categorical.
        """
        from ..io import write_zarr

        # TODO: What is a bool for chunks supposed to do?
        if isinstance(chunks, bool):
            msg = (
                "Passing `write_zarr(adata, chunks=True)` is no longer supported. "
                "Please pass `write_zarr(adata)` instead."
            )
            raise ValueError(msg)
        write_zarr(
            store,
            self,
            chunks=chunks,
            convert_strings_to_categoricals=convert_strings_to_categoricals,
        )

    def chunked_X(self, chunk_size: int | None = None):
        """\
        Return an iterator over the rows of the data matrix :attr:`X`.

        Parameters
        ----------
        chunk_size
            Row size of a single chunk.
        """
        if chunk_size is None:
            # Should be some adaptive code
            chunk_size = 6000
        start = 0
        n = self.n_obs
        for _ in range(int(n // chunk_size)):
            end = start + chunk_size
            yield (self.X[start:end], start, end)
            start = end
        if start < n:
            yield (self.X[start:n], start, n)

    @old_positionals("replace")
    def chunk_X(
        self,
        select: int | Sequence[int] | np.ndarray = 1000,
        *,
        replace: bool = True,
    ):
        """\
        Return a chunk of the data matrix :attr:`X` with random or specified indices.

        Parameters
        ----------
        select
            Depending on the type:

            :class:`int`
                A random chunk with `select` rows will be returned.
            :term:`sequence` (e.g. a list, tuple or numpy array) of :class:`int`
                A chunk with these indices will be returned.

        replace
            If `select` is an integer then `True` means random sampling of
            indices with replacement, `False` without replacement.
        """
        if isinstance(select, int):
            select = select if select < self.n_obs else self.n_obs
            choice = np.random.choice(self.n_obs, select, replace)
        elif isinstance(select, np.ndarray | Sequence):
            choice = np.asarray(select)
        else:
            msg = "select should be int or array"
            raise ValueError(msg)

        reverse = None
        if self.isbacked:
            # h5py can only slice with a sorted list of unique index values
            # so random batch with indices [2, 2, 5, 3, 8, 10, 8] will fail
            # this fixes the problem
            indices, reverse = np.unique(choice, return_inverse=True)
            selection = self.X[indices.tolist()]
        else:
            selection = self.X[choice]

        selection = selection.toarray() if issparse(selection) else selection
        return selection if reverse is None else selection[reverse]

    def _has_X(self) -> bool:
        """
        Check if X is None.

        This is more efficient than trying `adata.X is None` for views, since creating
        views (at least anndata's kind) can be expensive.
        """
        if not self.is_view:
            return self.X is not None
        else:
            return self._adata_ref.X is not None

    # --------------------------------------------------------------------------
    # all of the following is for backwards compat
    # --------------------------------------------------------------------------

    @property
    @deprecated("is_view")
    def isview(self):
        return self.is_view

    def _clean_up_old_format(self, uns):
        # multicolumn keys
        # all of the rest is only for backwards compat
        for bases in [["obs", "smp"], ["var"]]:
            axis = bases[0]
            for k in [f"{p}{base}_keys_multicol" for p in ["", "_"] for base in bases]:
                if uns and k in uns:
                    keys = list(uns[k])
                    del uns[k]
                    break
            else:
                keys = []
            # now, for compat, fill the old multicolumn entries into obsm and varm
            # and remove them from obs and var
            m_attr = getattr(self, f"_{axis}m")
            for key in keys:
                m_attr[key] = self._get_and_delete_multicol_field(axis, key)

    def _get_and_delete_multicol_field(self, a, key_multicol):
        keys = [k for k in getattr(self, a).columns if k.startswith(key_multicol)]
        values = getattr(self, a)[keys].values
        getattr(self, a).drop(keys, axis=1, inplace=True)
        return values


@AnnData._remove_unused_categories.register(Dataset2D)
@staticmethod
def _remove_unused_categories_xr(
    df_full: Dataset2D, df_sub: Dataset2D, uns: dict[str, Any]
):
    pass  # this is handled automatically by the categorical arrays themselves i.e., they dedup upon access.


def _check_2d_shape(X):
    """\
    Check shape of array or sparse matrix.

    Assure that X is always 2D: Unlike numpy we always deal with 2D arrays.
    """
    if X.dtype.names is None and len(X.shape) != 2:
        msg = f"X needs to be 2-dimensional, not {len(X.shape)}-dimensional."
        raise ValueError(msg)


def _infer_shape_for_axis(
    xxx: pd.DataFrame | Mapping[str, Iterable[Any]] | None,
    xxxm: np.ndarray | Mapping[str, Sequence[Any]] | None,
    layers: Mapping[str, np.ndarray | sparse.spmatrix] | None,
    xxxp: np.ndarray | Mapping[str, Sequence[Any]] | None,
    axis: Literal[0, 1],
) -> int | None:
    for elem in [xxx, xxxm, xxxp]:
        if elem is not None and hasattr(elem, "shape"):
            return elem.shape[0]
    for elem, id in zip([layers, xxxm, xxxp], ["layers", "xxxm", "xxxp"], strict=True):
        if elem is not None:
            elem = cast("Mapping", elem)
            for sub_elem in elem.values():
                if hasattr(sub_elem, "shape"):
                    size = cast("int", sub_elem.shape[axis if id == "layers" else 0])
                    return size
    return None


def _infer_shape(
    obs: pd.DataFrame | Mapping[str, Iterable[Any]] | None = None,
    var: pd.DataFrame | Mapping[str, Iterable[Any]] | None = None,
    *,
    obsm: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
    varm: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
    layers: Mapping[str, np.ndarray | sparse.spmatrix] | None = None,
    obsp: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
    varp: np.ndarray | Mapping[str, Sequence[Any]] | None = None,
):
    return (
        _infer_shape_for_axis(obs, obsm, layers, obsp, 0),
        _infer_shape_for_axis(var, varm, layers, varp, 1),
    )
