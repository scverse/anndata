from collections import OrderedDict, abc as cabc
from pathlib import Path
from typing import Any, MutableMapping, Union, List
from anndata._core.aligned_mapping import Layers, PairwiseArrays
from anndata._core.anndata import StorageType, _check_2d_shape
from anndata._core.index import Index
from anndata._core.raw import Raw
from anndata.compat import _move_adj_mtx
from anndata.utils import convert_to_dict

import zarr
import pandas as pd

from ..._core import AnnData, AxisArrays
from .utils import read_dispatched

class SingleDimensionAxisArraysRemote(AxisArrays):
    def __getitem__(self, key: str):
        return self._data[key][()]
    
    def __getattr__(self, __name: str):
        # If we a method has been accessed that is not here, try the pandas implementation
        if hasattr(pd.DataFrame, __name):
            return self.to_df().__getattribute__(__name)
        return object.__getattribute__(self, __name)



class AnnDataRemote(AnnData):

    def __init__(
        self,
        X = None,
        obs = None,
        var = None,
        uns = None,
        obsm = None,
        varm = None,
        layers = None,
        raw = None,
        dtype = None,
        shape = None,
        filename = None,
        filemode = None,
        asview = False,
        *,
        obsp,
        varp,
        oidx,
        vidx,
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

        # annotations
        self._obs = AxisArrays(self, 0, vals=convert_to_dict(obs))
        self._var = AxisArrays(self, 0, vals=convert_to_dict(var))

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
        self._obsm = AxisArrays(self, 0, vals=convert_to_dict(obsm))
        self._varm = AxisArrays(self, 1, vals=convert_to_dict(varm))

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


    def __eq__(self, other):
        """Equality testing"""
        raise NotImplementedError(
            "Equality comparisons are not supported for AnnData objects, "
            "instead compare the desired attributes."
        )

    @property
    def obs_names(self) -> pd.Index:
        """Names of observations (alias for `.obs.index`)."""
        return pd.Index(self.obs['_index'])

    @property
    def var_names(self) -> pd.Index:
        """Names of variables (alias for `.var.index`)."""
        return pd.Index(self.var['_index'])

    def obs_keys(self) -> List[str]:
        """List keys of observation annotation :attr:`obs`."""
        return self._obs.keys()

    def var_keys(self) -> List[str]:
        """List keys of variable annotation :attr:`var`."""
        return self._var.keys()

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

    # def obs_vector(self, k: str, *, layer: Optional[str] = None) -> np.ndarray:
    #     """\
    #     Convenience function for returning a 1 dimensional ndarray of values
    #     from :attr:`X`, :attr:`layers`\\ `[k]`, or :attr:`obs`.

    #     Made for convenience, not performance.
    #     Intentionally permissive about arguments, for easy iterative use.

    #     Params
    #     ------
    #     k
    #         Key to use. Should be in :attr:`var_names` or :attr:`obs`\\ `.columns`.
    #     layer
    #         What layer values should be returned from. If `None`, :attr:`X` is used.

    #     Returns
    #     -------
    #     A one dimensional nd array, with values for each obs in the same order
    #     as :attr:`obs_names`.
    #     """
    #     if layer == "X":
    #         if "X" in self.layers:
    #             pass
    #         else:
    #             warnings.warn(
    #                 "In a future version of AnnData, access to `.X` by passing"
    #                 " `layer='X'` will be removed. Instead pass `layer=None`.",
    #                 FutureWarning,
    #             )
    #             layer = None
    #     return get_vector(self, k, "obs", "var", layer=layer)

    # def var_vector(self, k, *, layer: Optional[str] = None) -> np.ndarray:
    #     """\
    #     Convenience function for returning a 1 dimensional ndarray of values
    #     from :attr:`X`, :attr:`layers`\\ `[k]`, or :attr:`obs`.

    #     Made for convenience, not performance. Intentionally permissive about
    #     arguments, for easy iterative use.

    #     Params
    #     ------
    #     k
    #         Key to use. Should be in :attr:`obs_names` or :attr:`var`\\ `.columns`.
    #     layer
    #         What layer values should be returned from. If `None`, :attr:`X` is used.

    #     Returns
    #     -------
    #     A one dimensional nd array, with values for each var in the same order
    #     as :attr:`var_names`.
    #     """
    #     if layer == "X":
    #         if "X" in self.layers:
    #             pass
    #         else:
    #             warnings.warn(
    #                 "In a future version of AnnData, access to `.X` by passing "
    #                 "`layer='X'` will be removed. Instead pass `layer=None`.",
    #                 FutureWarning,
    #             )
    #             layer = None
    # #     return get_vector(self, k, "var", "obs", layer=layer)

    # def to_memory(self, copy=True) -> "AnnData":
    #     """Return a new AnnData object with all backed arrays loaded into memory.

    #     Params
    #     ------
    #         copy:
    #             Whether the arrays that are already in-memory should be copied.

    #     Example
    #     -------

    #     .. code:: python

    #         import anndata
    #         backed = anndata.read_h5ad("file.h5ad", backed="r")
    #         mem = backed[backed.obs["cluster"] == "a", :].to_memory()
    #     """
    #     new = {}
    #     for attr_name in [
    #         "X",
    #         "obs",
    #         "var",
    #         "obsm",
    #         "varm",
    #         "obsp",
    #         "varp",
    #         "layers",
    #         "uns",
    #     ]:
    #         attr = getattr(self, attr_name, None)
    #         if attr is not None:
    #             new[attr_name] = to_memory(attr, copy)

    #     if self.raw is not None:
    #         new["raw"] = {
    #             "X": to_memory(self.raw.X, copy),
    #             "var": to_memory(self.raw.var, copy),
    #             "varm": to_memory(self.raw.varm, copy),
    #         }

    #     if self.isbacked:
    #         self.file.close()

    # #     return AnnData(**new)

    # def concatenate(
    #     self,
    #     *adatas: "AnnData",
    #     join: str = "inner",
    #     batch_key: str = "batch",
    #     batch_categories: Sequence[Any] = None,
    #     uns_merge: Optional[str] = None,
    #     index_unique: Optional[str] = "-",
    #     fill_value=None,
    # ) -> "AnnData":
    #     """\
    #     Concatenate along the observations axis.

    #     The :attr:`uns`, :attr:`varm` and :attr:`obsm` attributes are ignored.

    #     Currently, this works only in `'memory'` mode.

    #     .. note::

    #         For more flexible and efficient concatenation, see: :func:`~anndata.concat`.

    #     Parameters
    #     ----------
    #     adatas
    #         AnnData matrices to concatenate with. Each matrix is referred to as
    #         a “batch”.
    #     join
    #         Use intersection (`'inner'`) or union (`'outer'`) of variables.
    #     batch_key
    #         Add the batch annotation to :attr:`obs` using this key.
    #     batch_categories
    #         Use these as categories for the batch annotation. By default, use increasing numbers.
    #     uns_merge
    #         Strategy to use for merging entries of uns. These strategies are applied recusivley.
    #         Currently implemented strategies include:

    #         * `None`: The default. The concatenated object will just have an empty dict for `uns`.
    #         * `"same"`: Only entries which have the same value in all AnnData objects are kept.
    #         * `"unique"`: Only entries which have one unique value in all AnnData objects are kept.
    #         * `"first"`: The first non-missing value is used.
    #         * `"only"`: A value is included if only one of the AnnData objects has a value at this
    #           path.
    #     index_unique
    #         Make the index unique by joining the existing index names with the
    #         batch category, using `index_unique='-'`, for instance. Provide
    #         `None` to keep existing indices.
    #     fill_value
    #         Scalar value to fill newly missing values in arrays with. Note: only applies to arrays
    #         and sparse matrices (not dataframes) and will only be used if `join="outer"`.

    #         .. note::
    #             If not provided, the default value is `0` for sparse matrices and `np.nan`
    #             for numpy arrays. See the examples below for more information.

    #     Returns
    #     -------
    #     :class:`~anndata.AnnData`
    #         The concatenated :class:`~anndata.AnnData`, where `adata.obs[batch_key]`
    #         stores a categorical variable labeling the batch.

    #     Notes
    #     -----

    #     .. warning::

    #        If you use `join='outer'` this fills 0s for sparse data when
    #        variables are absent in a batch. Use this with care. Dense data is
    #        filled with `NaN`. See the examples.

    #     Examples
    #     --------
    #     Joining on intersection of variables.

    #     >>> adata1 = AnnData(
    #     ...     np.array([[1, 2, 3], [4, 5, 6]]),
    #     ...     dict(obs_names=['s1', 's2'], anno1=['c1', 'c2']),
    #     ...     dict(var_names=['a', 'b', 'c'], annoA=[0, 1, 2]),
    #     ... )
    #     >>> adata2 = AnnData(
    #     ...     np.array([[1, 2, 3], [4, 5, 6]]),
    #     ...     dict(obs_names=['s3', 's4'], anno1=['c3', 'c4']),
    #     ...     dict(var_names=['d', 'c', 'b'], annoA=[0, 1, 2]),
    #     ... )
    #     >>> adata3 = AnnData(
    #     ... np.array([[1, 2, 3], [4, 5, 6]]),
    #     ...     dict(obs_names=['s1', 's2'], anno2=['d3', 'd4']),
    #     ...     dict(var_names=['d', 'c', 'b'], annoA=[0, 2, 3], annoB=[0, 1, 2]),
    #     ... )
    #     >>> adata = adata1.concatenate(adata2, adata3)
    #     >>> adata
    #     AnnData object with n_obs × n_vars = 6 × 2
    #         obs: 'anno1', 'anno2', 'batch'
    #         var: 'annoA-0', 'annoA-1', 'annoA-2', 'annoB-2'
    #     >>> adata.X
    #     array([[2, 3],
    #            [5, 6],
    #            [3, 2],
    #            [6, 5],
    #            [3, 2],
    #            [6, 5]])
    #     >>> adata.obs
    #          anno1 anno2 batch
    #     s1-0    c1   NaN     0
    #     s2-0    c2   NaN     0
    #     s3-1    c3   NaN     1
    #     s4-1    c4   NaN     1
    #     s1-2   NaN    d3     2
    #     s2-2   NaN    d4     2
    #     >>> adata.var.T
    #              b  c
    #     annoA-0  1  2
    #     annoA-1  2  1
    #     annoA-2  3  2
    #     annoB-2  2  1

    #     Joining on the union of variables.

    #     >>> outer = adata1.concatenate(adata2, adata3, join='outer')
    #     >>> outer
    #     AnnData object with n_obs × n_vars = 6 × 4
    #         obs: 'anno1', 'anno2', 'batch'
    #         var: 'annoA-0', 'annoA-1', 'annoA-2', 'annoB-2'
    #     >>> outer.var.T
    #                a    b    c    d
    #     annoA-0  0.0  1.0  2.0  NaN
    #     annoA-1  NaN  2.0  1.0  0.0
    #     annoA-2  NaN  3.0  2.0  0.0
    #     annoB-2  NaN  2.0  1.0  0.0
    #     >>> outer.var_names
    #     Index(['a', 'b', 'c', 'd'], dtype='object')
    #     >>> outer.X
    #     array([[ 1.,  2.,  3., nan],
    #            [ 4.,  5.,  6., nan],
    #            [nan,  3.,  2.,  1.],
    #            [nan,  6.,  5.,  4.],
    #            [nan,  3.,  2.,  1.],
    #            [nan,  6.,  5.,  4.]])
    #     >>> outer.X.sum(axis=0)
    #     array([nan, 25., 23., nan])
    #     >>> import pandas as pd
    #     >>> Xdf = pd.DataFrame(outer.X, columns=outer.var_names)
    #     >>> Xdf
    #          a    b    c    d
    #     0  1.0  2.0  3.0  NaN
    #     1  4.0  5.0  6.0  NaN
    #     2  NaN  3.0  2.0  1.0
    #     3  NaN  6.0  5.0  4.0
    #     4  NaN  3.0  2.0  1.0
    #     5  NaN  6.0  5.0  4.0
    #     >>> Xdf.sum()
    #     a     5.0
    #     b    25.0
    #     c    23.0
    #     d    10.0
    #     dtype: float64

    #     One way to deal with missing values is to use masked arrays:

    #     >>> from numpy import ma
    #     >>> outer.X = ma.masked_invalid(outer.X)
    #     >>> outer.X
    #     masked_array(
    #       data=[[1.0, 2.0, 3.0, --],
    #             [4.0, 5.0, 6.0, --],
    #             [--, 3.0, 2.0, 1.0],
    #             [--, 6.0, 5.0, 4.0],
    #             [--, 3.0, 2.0, 1.0],
    #             [--, 6.0, 5.0, 4.0]],
    #       mask=[[False, False, False,  True],
    #             [False, False, False,  True],
    #             [ True, False, False, False],
    #             [ True, False, False, False],
    #             [ True, False, False, False],
    #             [ True, False, False, False]],
    #       fill_value=1e+20)
    #     >>> outer.X.sum(axis=0).data
    #     array([ 5., 25., 23., 10.])

    #     The masked array is not saved but has to be reinstantiated after saving.

    #     >>> outer.write('./test.h5ad')
    #     >>> from anndata import read_h5ad
    #     >>> outer = read_h5ad('./test.h5ad')
    #     >>> outer.X
    #     array([[ 1.,  2.,  3., nan],
    #            [ 4.,  5.,  6., nan],
    #            [nan,  3.,  2.,  1.],
    #            [nan,  6.,  5.,  4.],
    #            [nan,  3.,  2.,  1.],
    #            [nan,  6.,  5.,  4.]])

    #     For sparse data, everything behaves similarly,
    #     except that for `join='outer'`, zeros are added.

    #     >>> from scipy.sparse import csr_matrix
    #     >>> adata1 = AnnData(
    #     ...     csr_matrix([[0, 2, 3], [0, 5, 6]], dtype=np.float32),
    #     ...     dict(obs_names=['s1', 's2'], anno1=['c1', 'c2']),
    #     ...     dict(var_names=['a', 'b', 'c']),
    #     ... )
    #     >>> adata2 = AnnData(
    #     ...     csr_matrix([[0, 2, 3], [0, 5, 6]], dtype=np.float32),
    #     ...     dict(obs_names=['s3', 's4'], anno1=['c3', 'c4']),
    #     ...     dict(var_names=['d', 'c', 'b']),
    #     ... )
    #     >>> adata3 = AnnData(
    #     ... csr_matrix([[1, 2, 0], [0, 5, 6]], dtype=np.float32),
    #     ...     dict(obs_names=['s5', 's6'], anno2=['d3', 'd4']),
    #     ...     dict(var_names=['d', 'c', 'b']),
    #     ... )
    #     >>> adata = adata1.concatenate(adata2, adata3, join='outer')
    #     >>> adata.var_names
    #     Index(['a', 'b', 'c', 'd'], dtype='object')
    #     >>> adata.X.toarray()
    #     array([[0., 2., 3., 0.],
    #            [0., 5., 6., 0.],
    #            [0., 3., 2., 0.],
    #            [0., 6., 5., 0.],
    #            [0., 0., 2., 1.],
    #            [0., 6., 5., 0.]], dtype=float32)
    #     """
    #     from .merge import concat, merge_outer, merge_dataframes, merge_same

    #     warnings.warn(
    #         "The AnnData.concatenate method is deprecated in favour of the "
    #         "anndata.concat function. Please use anndata.concat instead.\n\n"
    #         "See the tutorial for concat at: "
    #         "https://anndata.readthedocs.io/en/latest/concatenation.html",
    #         FutureWarning,
    #     )

    #     if self.isbacked:
    #         raise ValueError("Currently, concatenate only works in memory mode.")

    #     if len(adatas) == 0:
    #         return self.copy()
    #     elif len(adatas) == 1 and not isinstance(adatas[0], AnnData):
    #         adatas = adatas[0]  # backwards compatibility
    #     all_adatas = (self,) + tuple(adatas)

    #     out = concat(
    #         all_adatas,
    #         axis=0,
    #         join=join,
    #         label=batch_key,
    #         keys=batch_categories,
    #         uns_merge=uns_merge,
    #         fill_value=fill_value,
    #         index_unique=index_unique,
    #         pairwise=False,
    #     )

    #     # Backwards compat (some of this could be more efficient)
    #     # obs used to always be an outer join
    #     out.obs = concat(
    #         [AnnData(sparse.csr_matrix(a.shape), obs=a.obs) for a in all_adatas],
    #         axis=0,
    #         join="outer",
    #         label=batch_key,
    #         keys=batch_categories,
    #         index_unique=index_unique,
    #     ).obs
    #     # Removing varm
    #     del out.varm
    #     # Implementing old-style merging of var
    #     if batch_categories is None:
    #         batch_categories = np.arange(len(all_adatas)).astype(str)
    #     pat = rf"-({'|'.join(batch_categories)})$"
    #     out.var = merge_dataframes(
    #         [a.var for a in all_adatas],
    #         out.var_names,
    #         partial(merge_outer, batch_keys=batch_categories, merge=merge_same),
    #     )
    #     out.var = out.var.iloc[
    #         :,
    #         (
    #             out.var.columns.str.extract(pat, expand=False)
    #             .fillna("")
    #             .argsort(kind="stable")
    #         ),
    #     ]

    #     return out

    # def var_names_make_unique(self, join: str = "-"):
    #     # Important to go through the setter so obsm dataframes are updated too
    #     self.var_names = utils.make_index_unique(self.var.index, join)

    # var_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    # def obs_names_make_unique(self, join: str = "-"):
    #     # Important to go through the setter so obsm dataframes are updated too
    #     self.obs_names = utils.make_index_unique(self.obs.index, join)

    # obs_names_make_unique.__doc__ = utils.make_index_unique.__doc__

    def __contains__(self, key: Any):
        raise AttributeError(
            "AnnData has no attribute __contains__, don;t check `in adata`."
        )

    # def _check_dimensions(self, key=None):
    #     if key is None:
    #         key = {"obs", "var", "obsm", "varm"}
    #     else:
    #         key = {key}
    #     if "obs" in key and len(self._obs) != self._n_obs:
    #         raise ValueError(
    #             "Observations annot. `obs` must have number of rows of `X`"
    #             f" ({self._n_obs}), but has {self._obs.shape[0]} rows."
    #         )
    #     if "var" in key and len(self._var) != self._n_vars:
    #         raise ValueError(
    #             "Variables annot. `var` must have number of columns of `X`"
    #             f" ({self._n_vars}), but has {self._var.shape[0]} rows."
    #         )
    #     if "obsm" in key:
    #         obsm = self._obsm
    #         if (
    #             not all([o.shape[0] == self._n_obs for o in obsm.values()])
    #             and len(obsm.dim_names) != self._n_obs
    #         ):
    #             raise ValueError(
    #                 "Observations annot. `obsm` must have number of rows of `X`"
    #                 f" ({self._n_obs}), but has {len(obsm)} rows."
    #             )
    #     if "varm" in key:
    #         varm = self._varm
    #         if (
    #             not all([v.shape[0] == self._n_vars for v in varm.values()])
    #             and len(varm.dim_names) != self._n_vars
    #         ):
    #             raise ValueError(
    #                 "Variables annot. `varm` must have number of columns of `X`"
    #                 f" ({self._n_vars}), but has {len(varm)} rows."
    #             )


def read_remote(store: Union[str, Path, MutableMapping, zarr.Group]) -> AnnData:
    if isinstance(store, Path):
        store = str(store)

    f = zarr.open(store, mode="r")

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith('/'):
            return AnnData(
                **{k: read_dispatched(v, callback) for k, v in elem.items()}
            )
        elif elem_name.startswith("raw."):
            return None
        elif elem_name in {"obs", "var"}:
            # override to only return AxisArray that will be accessed specially via our special AnnData object
            return {k: func(v) for k, v in elem.items()}
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata