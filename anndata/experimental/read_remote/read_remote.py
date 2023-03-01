from collections import OrderedDict, abc as cabc
from pathlib import Path
from typing import Any, MutableMapping, Union, List, Sequence
from anndata._core.aligned_mapping import Layers, PairwiseArrays
from anndata._core.anndata import StorageType, _check_2d_shape, _gen_dataframe
from anndata._core.file_backing import AnnDataFileManager
from anndata._core.index import Index
from anndata._core.raw import Raw
from anndata._io.specs.registry import read_elem
from anndata.compat import _move_adj_mtx, _read_attr
from anndata.utils import convert_to_dict

import zarr
import pandas as pd

from ..._core import AnnData, AxisArrays
from .utils import read_dispatched


# Initialization from something like
# CategoricalZarrArray(adata_local.obs['cell_type']['codes'].store, 'obs/cell_type/codes')
# Will need to work out the API better once I understand what the best practice here is.
class CategoricalZarrArray(zarr.core.Array):
    def __init__(self, group, *args, **kwargs):
        codes_path = group.path + "/codes"
        super().__init__(group.store.store, codes_path, *args, **kwargs)
        self.categories = group["categories"][()]
        self.ordered = bool(_read_attr(group.attrs, "ordered"))

    def __array__(self, *args):  # may need to override this, copied for now
        a = self[...]
        if args:
            a = a.astype(args[0])
        return a

    def __getitem__(self, selection):
        result = super().__getitem__(selection)
        return pd.Categorical.from_codes(
            codes=result,
            categories=self.categories,
            ordered=self.ordered,
        )


class AxisArraysRemote(AxisArrays):
    def __getattr__(self, __name: str):
        # If we a method has been accessed that is not here, try the pandas implementation
        if hasattr(pd.DataFrame, __name):
            return self.to_df().__getattribute__(__name)
        return object.__getattribute__(self, __name)

    def to_df(self) -> pd.DataFrame:
        """Convert to pandas dataframe."""
        df = pd.DataFrame(index=self.dim_names)
        for key in self.keys():
            z = self[key]
            if isinstance(z, zarr.Group) and "codes" in z:  # catrgoricql
                value = pd.Categorical.from_codes(
                    codes=read_elem(z["codes"]),
                    categories=read_elem(z["categories"]),
                    ordered=bool(_read_attr(z.attrs, "ordered")),
                )
            else:
                value = z[()]
            df[key] = value
        return df


class AnnDataRemote(AnnData):
    # TODO's here:
    # 1. Get an in-place copying system running
    # 2. Get a better sparse access pattern
    # 3. Re-write dataset with better chunking
    # 5. a `head` method

    def __init__(
        self,
        X=None,
        obs=None,
        var=None,
        uns=None,
        obsm=None,
        varm=None,
        layers=None,
        raw=None,
        dtype=None,
        shape=None,
        filename=None,
        filemode=None,
        asview=False,
        *,
        obsp=None,
        varp=None,
        oidx=None,
        vidx=None,
    ):
        # view attributes
        self._is_view = False
        self._adata_ref = None
        self._oidx = None
        self._vidx = None

        # ----------------------------------------------------------------------
        # various ways of initializing the data
        # ----------------------------------------------------------------------

        # check data type of X
        if filename is not None:
            self.file = AnnDataFileManager(self, filename, filemode)
        else:
            self.file = AnnDataFileManager(self, None)

        if X is not None:
            for s_type in StorageType:
                if isinstance(X, s_type.value):
                    break
            else:
                class_names = ", ".join(c.__name__ for c in StorageType.classes())
                raise ValueError(
                    f"`X` needs to be of one of {class_names}, not {type(X)}."
                )
            if shape is not None:
                raise ValueError("`shape` needs to be `None` if `X` is not `None`.")
            _check_2d_shape(X)
            # if type doesn’t match, a copy is made, otherwise, use a view
            if dtype is not None:
                X = X.astype(dtype)
            # data matrix and shape
            self._X = X
            self._n_obs, self._n_vars = self._X.shape
        else:
            self._X = None
            self._n_obs = len([] if obs is None else obs)
            self._n_vars = len([] if var is None else var)
            # check consistency with shape
            if shape is not None:
                if self._n_obs == 0:
                    self._n_obs = shape[0]
                else:
                    if self._n_obs != shape[0]:
                        raise ValueError("`shape` is inconsistent with `obs`")
                if self._n_vars == 0:
                    self._n_vars = shape[1]
                else:
                    if self._n_vars != shape[1]:
                        raise ValueError("`shape` is inconsistent with `var`")

        # annotations - need names already for AxisArrays to work.
        self.obs_names = pd.Index(
            (obs["index"] if "index" in obs else obs["_index"])[()]
        )
        self.var_names = pd.Index(
            (var["index"] if "index" in var else var["_index"])[()]
        )
        self._obs = AxisArraysRemote(self, 0, vals=convert_to_dict(obs))
        self._var = AxisArraysRemote(self, 1, vals=convert_to_dict(var))

        # now we can verify if indices match!
        # for attr_name, x_name, idx in x_indices:
        #     attr = getattr(self, attr_name)
        #     if isinstance(attr.index, pd.RangeIndex):
        #         attr.index = idx
        #     elif not idx.equals(attr.index):
        #         raise ValueError(f"Index of {attr_name} must match {x_name} of X.")

        # unstructured annotations
        self.uns = uns or OrderedDict()

        # TODO: Think about consequences of making obsm a group in hdf
        self._obsm = AxisArraysRemote(self, 0, vals=convert_to_dict(obsm))
        self._varm = AxisArraysRemote(self, 1, vals=convert_to_dict(varm))

        self._obsp = PairwiseArrays(self, 0, vals=convert_to_dict(obsp))
        self._varp = PairwiseArrays(self, 1, vals=convert_to_dict(varp))

        # Backwards compat for connectivities matrices in uns["neighbors"]
        _move_adj_mtx({"uns": self._uns, "obsp": self._obsp})

        # self._check_dimensions()
        # self._check_uniqueness()

        if self.filename:
            assert not isinstance(
                raw, Raw
            ), "got raw from other adata but also filename?"
            if {"raw", "raw.X"} & set(self.file):
                raw = dict(X=None, **raw)
        if not raw:
            self._raw = None
        elif isinstance(raw, cabc.Mapping):
            self._raw = Raw(self, **raw)
        else:  # is a Raw from another AnnData
            self._raw = Raw(self, raw._X, raw.var, raw.varm)

        # clean up old formats
        self._clean_up_old_format(uns)

        # layers
        self._layers = Layers(self, layers)

    # TODO: this is not quite complete in the original but also here, what do we do about this?
    def __delitem__(self, index: Index):
        obs, var = self._normalize_indices(index)
        # TODO: does this really work?
        if not self.isbacked:
            del self._X[obs, var]
        else:
            X = self.file["X"]
            del X[obs, var]
            self._set_backed("X", X)

    @property
    def obs_names(self) -> pd.Index:
        return self._obs_names

    @property
    def var_names(self) -> pd.Index:
        return self._var_names

    @obs_names.setter
    def obs_names(self, names: Sequence[str]):
        self._obs_names = names

    @var_names.setter
    def var_names(self, names: Sequence[str]):
        self._var_names = names


def read_remote(store: Union[str, Path, MutableMapping, zarr.Group]) -> AnnData:
    if isinstance(store, Path):
        store = str(store)

    f = zarr.open_consolidated(store, mode="r")

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            return AnnDataRemote(
                **{k: read_dispatched(v, callback) for k, v in elem.items()}
            )
        elif elem_name.startswith("/raw"):
            return None
        elif elem_name in {"/obs", "/var"}:
            # override to only return AxisArray that will be accessed specially via our special AnnData object
            return {k: read_dispatched(v, callback) for k, v in elem.items()}
        elif iospec.encoding_type == "categorical":
            return CategoricalZarrArray(elem)
        elif iospec.encoding_type in {"array", "string_array"}:
            return elem
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
