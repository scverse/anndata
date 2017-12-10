"""Main class and helper functions.
"""
import os
import h5py
from collections import Mapping, Sequence
from collections import OrderedDict
from enum import Enum
import numpy as np
from numpy import ma
import pandas as pd
from pandas.core.index import RangeIndex
from scipy import sparse as sp
from scipy.sparse.sputils import IndexMixin
from textwrap import dedent


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sp.spmatrix

    @classmethod
    def classes(cls):
        return tuple(c.value for c in cls.__members__.values())


class BoundRecArr(np.recarray):
    """A `np.recarray` to which fields can be added using `.['key']`.

    To enable this, it is bound to a instance of AnnData.

    Parameters
    ----------
    source : `np.ndarray`
        A (structured) numpy array.
    parent : object
        Any object to which the BoundRecArr shall be bound to.
    attr : string
        The name of the attribute as which it appears in parent.
    """
    _attr_choices = ['obsm', 'varm']

    def __new__(cls, input_array, parent, attr):
        arr = np.asarray(input_array).view(cls)
        arr._parent = parent
        arr._attr = attr
        return arr

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self._parent = getattr(obj, '_parent', None)
        self._attr = getattr(obj, '_attr', None)

    def copy(self):
        new = super(BoundRecArr, self).copy()
        new._parent = self._parent
        return new

    def flipped(self):
        new_attr = (self._attr_choices[1] if self._attr == self._attr_choices[0]
                    else self._attr_choices[0])
        return BoundRecArr(self, self._parent, new_attr)

    def keys(self):
        return self.dtype.names

    def __setitem__(self, key, arr):
        if arr.ndim == 1:
            raise ValueError('Use adata.obs or adata.var for 1-dimensional arrays.')
        if self.shape[0] != arr.shape[0]:
            raise ValueError('Can only assign an array of same length ({}), '
                             'not of length {}.'
                             .format(self.shape[0], arr.shape[0]))
        # the following always allocates a new array
        # even if the key already exists and dimensions match
        # TODO: one could check for this case
        # dtype
        merged_dtype = []
        found_key = False
        for descr in self.dtype.descr:
            if descr[0] == key:
                merged_dtype.append((key, arr.dtype, arr.shape[1]))
                found_key = True
            else:
                merged_dtype.append(descr)
        if not found_key:
            merged_dtype.append((key, arr.dtype, arr.shape[1]))
        # create new array
        new = np.empty(len(self), dtype=merged_dtype)
        # fill the array
        for name in new.dtype.names:
            if name == key:
                new[name] = arr
            else:
                new[name] = self[name]
        # make it a BoundRecArr
        # TODO: why can we not do this step before filling the array?
        new = BoundRecArr(new, self._parent, self._attr)
        setattr(self._parent, self._attr, new)

    def to_df(self):
        """Convert to pandas dataframe."""
        df = pd.DataFrame(index=RangeIndex(0, self.shape[0], name=None))
        for key in self.keys():
            value = self[key]
            for icolumn, column in enumerate(value.T):
                df['{}{}'.format(key, icolumn+1)] = column
        return df


# for backwards compat
def _find_corresponding_multicol_key(key, keys_multicol):
    """Find the corresponding multicolumn key."""
    for mk in keys_multicol:
        if key.startswith(mk) and 'of' in key:
            return mk
    return None


# for backwards compat
def _gen_keys_from_multicol_key(key_multicol, n_keys):
    """Generates single-column keys from multicolumn key."""
    keys = [('{}{:03}of{:03}')
            .format(key_multicol, i+1, n_keys) for i in range(n_keys)]
    return keys


def df_to_records_fixed_width(df):
    from pandas.api.types import is_string_dtype, is_categorical
    uns = {}  # unstructured dictionary for storing categories
    names = ['index']
    if is_string_dtype(df.index):
        max_len_index = df.index.map(len).max()
        index = df.index.values.astype('S{}'.format(max_len_index))
    else:
        index = df.index.values
    arrays = [index]
    for k in df.columns:
        names.append(k)
        if is_string_dtype(df[k]):
            c = df[k].astype('category')
            # transform to category only if there are less categories than
            # observations
            if len(c.cat.categories) < len(c):
                uns[k + '_categories'] = c.cat.categories.values
                arrays.append(c.cat.codes)
            else:
                max_len_index = df[k].map(len).max()
                arrays.append(df[k].values.astype('S{}'.format(max_len_index)))
        elif is_categorical(df[k]):
            uns[k + '_categories'] = df[k].cat.categories
            arrays.append(df[k].cat.codes)
        else:
            arrays.append(df[k].values)
    formats = [v.dtype for v in arrays]
    return np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats}), uns


def _fix_shapes(X, single_col=False):
    """Fix shapes of array or sparse matrix.

    Flatten 2d array with a single row or column to emulate numpys
    behavior of slicing.
    """
    if X.dtype.names is None and len(X.shape) not in {0, 1, 2}:
        raise ValueError('X needs to be 2-dimensional, not '
                         '{}-dimensional.'.format(len(X.shape)))
    if len(X.shape) == 2:
        # TODO: int immutable, copy of references to ints in X.shape
        # only valid until accidental change
        n_obs, n_vars = X.shape
        if n_obs == 1 and n_vars == 1:
            X = X[0, 0]
        elif n_obs == 1 or n_vars == 1:
            if sp.issparse(X): X = X.toarray()
            X = X.flatten()
    elif len(X.shape) == 1:
        n_vars = X.shape[0]
        n_obs = 1
    elif len(X.shape) == 1 and single_col:
        n_obs = X.shape[0]
        n_vars = 1
    else:
        n_vars = 1
        n_obs = 1
    return X, n_obs, n_vars


def _gen_dataframe(anno, length, index_names):
    if anno is None or len(anno) == 0:
        _anno = pd.DataFrame(index=RangeIndex(0, length, name=None).astype(str))
    else:
        for index_name in index_names:
            if index_name in anno:
                _anno = pd.DataFrame(
                    anno, index=anno[index_name],
                    columns=[k for k in anno.keys() if k != index_name])
                break
        else:
            _anno = pd.DataFrame(anno)
    return _anno


def save_sparse_csr(X, key='X'):
    from scipy.sparse.csr import csr_matrix
    X = csr_matrix(X)
    key_csr = key + '_csr'
    return {key_csr + '_data': X.data,
            key_csr + '_indices': X.indices,
            key_csr + '_indptr': X.indptr,
            key_csr + '_shape': np.array(X.shape)}


def load_sparse_csr(d, key='X'):
    from scipy.sparse.csr import csr_matrix
    key_csr = key + '_csr'
    d[key] = csr_matrix((d[key_csr + '_data'],
                         d[key_csr + '_indices'],
                         d[key_csr + '_indptr']),
                        shape=d[key_csr + '_shape'])
    del d[key_csr + '_data']
    del d[key_csr + '_indices']
    del d[key_csr + '_indptr']
    del d[key_csr + '_shape']
    return d


def read_h5ad(filename, init_filename=True):
    """Read `.h5ad`-formatted hdf5 file.

    Parameters
    ----------
    filename : `str`
        File name of data file.
    init_filename : `bool`, optional (default: `True`)
        If `True`, `filename` initializes the ``.filename`` attribute of the
        returned object and the data is *not* loaded into memory.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """
    if init_filename:
        return AnnData(filename=filename)
    else:
        d = _read_h5ad(filename)
        return AnnData(d)


def _read_h5ad(filename, attr=None, skip_attrs=None):
    """Return a dict with arrays for initializing AnnData.

    Parameters
    ----------
    filename : `str`
    attr : `str`
    skip_attrs : `tuple` of `str`
    """
    def postprocess_reading(key, value):
        if value.ndim == 1 and len(value) == 1:
            value = value[0]
        if value.dtype.kind == 'S':
            value = value.astype(str)
            # recover a dictionary that has been stored as a string
            if len(value) > 0:
                if value[0] == '{' and value[-1] == '}': value = eval(value)
        if (key != 'obs' and key != 'var' and key != '_obs' and key != '_var'
            and not isinstance(value, dict) and value.dtype.names is not None):
            # TODO: come up with a better way of solving this, see also below
            new_dtype = [((dt[0], 'U{}'.format(int(int(dt[1][2:])/4)))
                          if dt[1][1] == 'S' else dt) for dt in value.dtype.descr]
            value = value.astype(new_dtype)
        return key, value

    filename = str(filename)  # allow passing pathlib.Path objects
    d = {}
    with h5py.File(filename, 'r') as f:
        for key in f.keys():
            if attr is not None:
                if key != attr and not key.startswith(attr + '_csr'):
                    continue
            if skip_attrs is not None:
                if key in skip_attrs or key.startswith(skip_attrs):
                    continue
            # the '()' means 'read everything' (by contrast, ':' only works
            # if not reading a scalar type)
            value = f[key][()]
            key, value = postprocess_reading(key, value)
            d[key] = value
    csr_keys = [key.replace('_csr_data', '')
                for key in d if '_csr_data' in key]
    for key in csr_keys: d = load_sparse_csr(d, key=key)
    if attr is not None: return d[attr]

    return d


def write_h5ad(filename, adata,
                  compression='gzip', compression_opts=None):
    _write_h5ad(filename, adata=adata,
                   compression=compression, compression_opts=compression_opts)

# attr_value needs to be a dict
def _write_h5ad(filename, adata=None, attr_value=None,
                   compression='gzip', compression_opts=None):
    def preprocess_writing(value):
        if isinstance(value, dict):
            # hack for storing dicts
            value = np.array([str(value)])
        else:
            value = np.array(value)
            if value.ndim == 0: value = np.array([value])
        # make sure string format is chosen correctly
        if value.dtype.kind == 'U': value = value.astype(np.string_)
        return value
    filename = str(filename)  # allow passing pathlib.Path objects
    if not filename.endswith(('.h5', '.h5ad')):
        raise ValueError('Filename needs to end with \'.h5ad\'.')
    if not os.path.exists(os.path.dirname(filename)):
        # logg.info('creating directory', directory + '/', 'for saving output files')
        os.makedirs(directory)
    if attr_value is None:
        d = adata._to_dict_fixed_width_arrays()
    else:
        d = attr_value
    d_write = {}
    for key, value in d.items():
        if sp.issparse(value):
            for k, v in save_sparse_csr(value, key=key).items():
                d_write[k] = v
        else:
            d_write[key] = preprocess_writing(value)
            # some output about the data to write
            # print(type(value), value.dtype, value.dtype.kind, value.shape)
    with h5py.File(filename, 'w' if attr_value is None else 'r+') as f:
        for key, value in d_write.items():
            try:
                # ignore arrays with empty dtypes
                if value.dtype.descr:
                    f.create_dataset(key, data=value,
                                     compression=compression,
                                     compression_opts=compression_opts)
            except TypeError:
                # try writing it as byte strings
                try:
                    if value.dtype.names is None:
                        f.create_dataset(key, data=value.astype('S'),
                                         compression=compression,
                                         compression_opts=compression_opts)
                    else:
                        new_dtype = [(dt[0], 'S{}'.format(int(dt[1][2:])*4))
                                     for dt in value.dtype.descr]
                        f.create_dataset(key, data=value.astype(new_dtype),
                                         compression=compression,
                                         compression_opts=compression_opts)
                except Exception as e:
                    # logg.info(str(e))
                    warnings.warn('Could not save field with key = "{}" to hdf5 file.'
                                  .format(key))


class AnnData(IndexMixin):

    _main_narrative = dedent("""\
        :class:`~anndata.AnnData` stores a data matrix ``.X`` together with
        annotations of observations ``.obs``, variables ``.var`` and
        unstructured annotations ``.uns``. If setting an `.h5ad` backing file
        ``.filename``, data remains on the disk but is automatically loaded into
        memory if needed.

        .. raw:: html

            <img src="http://falexwolf.de/img/scanpy/anndata.svg"
                 style="width: 350px; margin: 10px 0px 15px 20px">

        An :class:`~anndata.AnnData` object ``adata`` can be sliced like a
        pandas dataframe, for example, ``adata_subset = adata[:,
        list_of_variable_names]``. Observation and variable names can be
        accessed via ``.obs_names`` and ``.var_names``, respectively.
        :class:`~anndata.AnnData`'s basic structure is similar to R's
        ExpressionSet [Huber15]_.

        If you find the `anndata` package useful, please consider citing [Wolf17]_.
        """)

    _example_concatenate = dedent("""\
        >>> adata1 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'obs_names': ['o1', 'o2'],
        >>>                   'anno1': ['c1', 'c2']},
        >>>                  {'var_names': ['a', 'b', 'c']})
        >>> adata2 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'obs_names': ['o3', 'o4'],
        >>>                   'anno1': ['c3', 'c4']},
        >>>                  {'var_names': ['b', 'c', 'd']})
        >>> adata3 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
        >>>                  {'obs_names': ['o5', 'o6'],
        >>>                   'anno2': ['d3', 'd4']},
        >>>                  {'var_names': ['b', 'c', 'd']})
        >>>
        >>> adata = adata1.concatenate([adata2, adata3])
        >>> adata.X
        [[ 2.  3.]
         [ 5.  6.]
         [ 1.  2.]
         [ 4.  5.]
         [ 1.  2.]
         [ 4.  5.]]
        >>> adata.obs
           anno1 anno2 batch
        o1    c1   NaN     0
        o2    c2   NaN     0
        o3    c3   NaN     1
        o4    c4   NaN     1
        o5   NaN    d3     2
        o6   NaN    d4     2
        """)

    def __init__(self, X=None, obs=None, var=None, uns=None,
                 obsm=None, varm=None,
                 filename=None,
                 dtype='float32', single_col=False):

        # init backing filename (might be `None`)
        self._filename = filename

        # init from file
        if filename is not None:
            if any((X, obs, var, uns, obsm, varm)):
                raise ValueError(
                    'If initializing from `filename`, '
                    'no further arguments may be passed.')
            # returns dictionary
            X = _read_h5ad(filename, skip_attrs=('X',))

        # generate from dictionary
        if isinstance(X, Mapping):
            if any((obs, var, uns, obsm, varm)):
                raise ValueError(
                    'If `X` is a dict no further arguments must be provided.')
            X, obs, var, uns, obsm, varm = self._from_dict(X)

        # check data type of X
        if X is not None:
            for s_type in StorageType:
                if isinstance(X, s_type.value):
                    break
            else:
                class_names = ', '.join(c.__name__ for c in StorageType.classes())
                raise ValueError('`X` needs to be of one of {}, not {}.'
                                 .format(class_names, type(X)))
            # if type doesn't match, a copy is made, otherwise, use a view
            if sp.issparse(X) or isinstance(X, ma.MaskedArray):
                # TODO: maybe use view on data attribute of sparse matrix
                #       as in readwrite.read_10x_h5
                if X.dtype != np.dtype(dtype): X = X.astype(dtype)
            else:  # is np.ndarray
                X = X.astype(dtype, copy=False)
            # data matrix and shape
            self._X, self._n_obs, self._n_vars = _fix_shapes(X)
        else:
            self._n_obs = len(obs)
            self._n_vars = len(var)

        # annotations
        self._obs = _gen_dataframe(obs, self._n_obs,
                                   ['obs_names', 'row_names', 'smp_names'])
        self._var = _gen_dataframe(var, self._n_vars, ['var_names', 'col_names'])

        # unstructured annotations
        self._uns = uns or {}

        # multi-dimensional array annotations
        if obsm is None: obsm = np.empty(self._n_obs, dtype=[])
        if varm is None: varm = np.empty(self._n_vars, dtype=[])
        self._obsm = BoundRecArr(obsm, self, 'obsm')
        self._varm = BoundRecArr(varm, self, 'varm')

        self._check_dimensions()

        # clean up old formats
        self._clean_up_old_format(uns)

    def __repr__(self):
        if self.filename is not None:
            backed_at = 'backed at \'{}\''.format(self.filename)
        else:
            backed_at = ''
        descr = (
            'AnnData object with n_obs × n_vars = {} × {} {}\n'
            '    obs_keys = {}\n'
            '    var_keys = {}\n'
            '    uns_keys = {}\n'
            '    obsm_keys = {}\n'
            '    varm_keys = {}'
            .format(self._n_obs, self._n_vars, backed_at,
                    self.obs_keys(), self.var_keys(),
                    self.uns_keys(),
                    self.obsm_keys(), self.varm_keys()))
        return descr

    def _read(self, attr):
        return _read_h5ad(self.filename, attr)

    def _write(self, attr, value):
        return _write_h5ad(self.filename, attr_value={attr: value})

    @property
    def filename(self):
        """Filename for backing object as `.h5ad` file on disk."""
        return self._filename

    @filename.setter
    def filename(self, value):
        write_h5ad(value, self)
        self._filename = value
        self._X = None

    @property
    def shape(self):
        """Shape of data matrix: (n_obs, n_vars)."""
        return self.n_obs, self.n_vars

    # for backwards compat
    @property
    def data(self):
        """Deprecated access to X."""
        print('DEPRECATION WARNING: use attribute `.X` instead of `.data`, '
              '`.data` will be removed in the future.')
        return self._X

    # for backwards compat
    @data.setter
    def data(self, value):
        print('DEPRECATION WARNING: use attribute `.X` instead of `.data`, '
              '`.data` will be removed in the future.')
        if not self.shape == value.shape:
            raise ValueError('Data matrix has wrong shape {}, need to be {}'
                             .format(value.shape, self.shape))
        self._X = value

    @property
    def X(self):
        """Data matrix of shape `n_obs` × `n_vars` (`np.ndarray`, `sp.spmatrix`)."""
        if self.filename is None: return self._X
        else: return self._read('X')

    @X.setter
    def X(self, value):
        if not self.shape == value.shape:
            raise ValueError('Data matrix has wrong shape {}, need to be {}'
                             .format(value.shape, self.shape))
        if self.filename is None: self._X = value
        else: self._write('X', value)

    @property
    def n_smps(self):
        """Deprecated: Number of observations."""
        return self._n_obs

    @property
    def n_obs(self):
        """Number of observations (rows)."""
        return self._n_obs

    @property
    def n_vars(self):
        """Number of variables / features."""
        return self._n_vars

    @property
    def obs(self):
        """One-dimensional annotation of observations (`pd.DataFrame`)."""
        return self._obs

    # for backwards compat
    @property
    def smp(self):
        """Deprecated: One-dimensional annotation of observations (`pd.DataFrame`)."""
        return self._obs

    @obs.setter
    def obs(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        self._obs = value

    # for backwards compat
    @smp.setter
    def smp(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        self._obs = value

    @property
    def var(self):
        """One-dimensional annotation of variables/ features (`pd.DataFrame`)."""
        return self._var

    @var.setter
    def var(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_vars:
            raise ValueError('Length does not match.')
        self._var = value

    @property
    def uns(self):
        """Unstructured annotation (ordered dictionary)."""
        return self._uns

    @uns.setter
    def uns(self, value):
        self._uns = value

    @property
    def add(self):
        """Deprecated, remains for backwards compatibility."""
        # FutureWarning and DeprecationWarning are not visible by default, use print
        print('DEPRECATION WARNING: use attribute `.uns` instead of `.add`, '
              '`.add` will be removed in the future.')
        return self._uns

    @add.setter
    def add(self, value):
        print('DEPRECATION WARNING: use attribute `.uns` instead of `.add`, '
              '`.add` will be removed in the future.')
        self._uns = value

    @property
    def obsm(self):
        """Multi-dimensional annotation of observations (mutable structured `np.ndarray`).

        Stores for each key, a two or higher-dimensional `np.ndarray` of length
        `n_obs`. Is sliced with `data` and `obs` but behaves otherwise like a
        `dict`.
        """
        return self._obsm

    # for backwards compat
    @property
    def smpm(self):
        """Multi-dimensional annotation of observations (mutable structured `np.ndarray`).

        Stores for each key, a two or higher-dimensional `np.ndarray` of length
        `n_obs`. Is sliced with `data` and `obs` but behaves otherwise like a
        `dict`.
        """
        return self._obsm

    @obsm.setter
    def obsm(self, value):
        if not isinstance(value, np.ndarray):
            raise ValueError('Can only assign np.ndarray.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        value = BoundRecArr(value, self, 'obsm')
        self._obsm = value

    # for backwards compat
    @smpm.setter
    def smpm(self, value):
        if not isinstance(value, np.ndarray):
            raise ValueError('Can only assign np.ndarray.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        value = BoundRecArr(value, self, 'obsm')
        self._obsm = value

    @property
    def varm(self):
        """Multi-dimensional annotation of variables/ features (mutable structured `np.ndarray`).

        Stores for each key, a two or higher-dimensional `np.ndarray` of length
        `n_vars`. Is sliced with `data` and `var` but behaves otherwise like a
        `dict`.
        """
        return self._varm

    @varm.setter
    def varm(self, value):
        if not isinstance(value, np.ndarray):
            raise ValueError('Can only assign np.ndarray.')
        if len(value) != self.n_vars:
            raise ValueError('Length does not match.')
        value = BoundRecArr(value, self, 'varm')
        self._varm = value

    @property
    def obs_names(self):
        """Index for observations (`obs.index`)."""
        return self._obs.index.values

    # backwards compat
    @property
    def smp_names(self):
        """Index for observations (`obs.index`)."""
        return self._obs.index.values

    @obs_names.setter
    def obs_names(self, names):
        self._obs.index = names

    # backwards compat
    @smp_names.setter
    def smp_names(self, names):
        self._obs.index = names

    @property
    def var_names(self):
        """Index for variables (`var.index`)."""
        return self._var.index.values

    @var_names.setter
    def var_names(self, names):
        self._var.index = names

    def obs_keys(self):
        """List keys of observation annotation `obs`."""
        return self._obs.keys().tolist()

    # for backwards compat
    def smp_keys(self):
        """List keys of observation annotation `obs`."""
        return self._obs.keys().tolist()

    def var_keys(self):
        """List keys of variable annotation `var`."""
        return self._var.keys().tolist()

    def obsm_keys(self):
        """List keys of observation annotation `obsm`."""
        return list(self._obsm.keys())

    # for backwards compat
    def smpm_keys(self):
        """List keys of observation annotation `obsm`."""
        return list(self._obsm.keys())

    def varm_keys(self):
        """List keys of variable annotation `varm`."""
        return list(self._varm.keys())

    def uns_keys(self):
        """List keys of unstructured annotation."""
        return sorted(list(self._uns.keys()))

    def _normalize_indices(self, packed_index):
        # deal with slicing with pd.Series
        if isinstance(packed_index, pd.Series):
            packed_index = packed_index.values
        if isinstance(packed_index, tuple) and len(packed_index) == 2:
            if isinstance(packed_index[1], pd.Series):
                packed_index = packed_index[0], packed_index[1].values
            if isinstance(packed_index[0], pd.Series):
                packed_index = packed_index[0].values, packed_index[1]
        obs, var = super(AnnData, self)._unpack_index(packed_index)
        obs = self._normalize_index(obs, self.obs_names)
        var = self._normalize_index(var, self.var_names)
        return obs, var

    def _normalize_index(self, index, names):
        def name_idx(i):
            if isinstance(i, str):
                # `where` returns an 1-tuple (1D array) of found indices
                i = np.where(names == i)[0]
                if len(i) == 0:  # returns array of length 0 if nothing is found
                    raise IndexError('Name "{}" is not valid variable or observation index.'
                                     .format(index))
                i = i[0]
            return i

        if isinstance(index, slice):
            start = name_idx(index.start)
            stop = name_idx(index.stop)
            # string slices can only be inclusive, so +1 in that case
            if isinstance(index.stop, str):
                stop = None if stop is None else stop + 1
            step = index.step
        elif isinstance(index, (int, str)):
            start = name_idx(index)
            stop = start + 1
            step = 1
        elif isinstance(index, (Sequence, np.ndarray)):
            return np.fromiter(map(name_idx, index), 'int64')
        else:
            raise IndexError('Unknown index {!r} of type {}'
                             .format(index, type(index)))

        return slice(start, stop, step)

    # TODO: this is not quite complete...
    def __delitem__(self, index):
        obs, var = self._normalize_indices(index)
        # TODO: does this really work?
        if self.filename is None:
            del self._X[obs, var]
        else:
            X = self._read('X')
            del X[obs, var]
            self._write('X', X)
        if var == slice(None):
            del self._obs.iloc[obs, :]
        if obs == slice(None):
            del self._var.iloc[var, :]

    def __getitem__(self, index):
        # Note: this cannot be made inplace
        # http://stackoverflow.com/questions/31916617/using-keyword-arguments-in-getitem-method-in-python
        obs, var = self._normalize_indices(index)
        if self.filename is None: X = self._X[obs, var]
        else: X = self._read('X')
        obs_new = self._obs.iloc[obs]
        obsm_new = self._obsm[obs]
        var_new = self._var.iloc[var]
        varm_new = self._varm[var]
        assert obs_new.shape[0] == X.shape[0], (obs, obs_new)
        assert var_new.shape[0] == X.shape[1], (var, var_new)
        uns_new = self._uns.copy()
        # slice sparse spatrices of n_obs × n_obs in self._uns
        if not (isinstance(obs, slice) and
                obs.start is None and obs.step is None and obs.stop is None):
            raised_warning = False
            for k, v in self._uns.items():
                if isinstance(v, sp.spmatrix) and v.shape == (self._n_obs, self._n_obs):
                    uns_new[k] = v.tocsc()[:, obs].tocsr()[obs, :]
        return AnnData(X, obs_new, var_new, uns_new, obsm_new, varm_new)

    def _inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[:, index], but inplace.
        """
        if self.filename is None:
            self._X = self._X[:, index]
            self._n_vars = self._X.shape[1]
        else:
            X = self._read('X')
            X = X[:, index]
            self._n_vars = X.shape[1]
            self._write('X', X)
        self._var = self._var.iloc[index]
        # TODO: the following should not be necessary!
        self._varm = BoundRecArr(self._varm[index], self, 'varm')
        return None

    def _inplace_subset_obs(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[index, :], but inplace.
        """
        if self.filename is None:
            self._X = self._X[index, :]
            self._n_obs = self._X.shape[0]
        else:
            X = self._read('X')
            X = X[index, :]
            self._n_obs = X.shape[0]
            self._write('X', X)
        for k, v in self._uns.items():
            if isinstance(v, sp.spmatrix) and v.shape == (self._n_obs, self._n_obs):
                self._uns[k] = v.tocsc()[:, index].tocsr()[index, :]
        self._obs = self._obs.iloc[index]
        # TODO: the following should not be necessary!
        self._obsm = BoundRecArr(self._obsm[index], self, 'obsm')
        return None

    def _get_obs_array(self, k):
        """Get an array along the observation dimension by first looking up
        obs_keys and then var_names."""
        x = (self._obs[k] if k in self.obs_keys()
             else self[:, k].data if k in set(self.var_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in obs_keys or var_names.'
                             .format(k))
        return x

    def _get_var_array(self, k):
        """Get an array along the variables dimension by first looking up
        var_keys and then obs_names."""
        x = (self._var[k] if k in self.var_keys()
             else self[k] if k in set(self.obs_names)
             else None)
        if x is None:
            raise ValueError('Did not find {} in var_keys or obs_names.'
                             .format(k))
        return x

    def __setitem__(self, index, val):
        obs, var = self._normalize_indices(index)
        if self.filename is None:
            self._X[obs, var] = val
        else:
            X = self._read('X')
            X[obs, var] = val
            self._write('X', X)

    def __len__(self):
        return self.shape[0]

    def transpose(self):
        """Transpose whole object.

        Data matrix is transposed, observations and variables are interchanged.
        """
        if self.filename is None: X = self._X
        else: X = _read('X')
        if sp.isspmatrix_csr(X):
            return AnnData(X.T.tocsr(), self._var, self._obs, self._uns,
                           self._varm.flipped(), self._obsm.flipped(),
                           filename=self.filename)
        return AnnData(X.T, self._var, self._obs, self._uns,
                       self._varm.flipped(), self._obsm.flipped(),
                       filename=self.filename)

    T = property(transpose)

    def copy(self):
        """Full copy with memory allocated."""
        if self.filename is None: X = self._X
        else: X = _read('X')
        return AnnData(X.copy() if self.filename is None else X,
                       self._obs.copy(),
                       self._var.copy(), self._uns.copy(),
                       self._obsm.copy(), self._varm.copy())

    def concatenate(self, adatas, batch_key='batch', batch_categories=None):
        if isinstance(adatas, AnnData): adatas = [adatas]
        joint_variables = self.var_names
        for adata2 in adatas:
            joint_variables = np.intersect1d(
                joint_variables, adata2.var_names, assume_unique=True)
        adatas_to_concat = []
        if batch_categories is None:
            categories = [str(i) for i in range(len(adatas)+1)]
        elif len(batch_categories) == len(adatas)+1:
            categories = batch_categories
        else:
            raise ValueError('Provide as many `batch_categories` as `adatas`.')
        for i, ad in enumerate([self] + adatas):
            ad = ad[:, joint_variables]
            ad.obs[batch_key] = pd.Categorical(
                ad.n_obs*[categories[i]], categories=categories)
            adatas_to_concat.append(ad)
        Xs = [ad.X for ad in adatas_to_concat]
        if sp.issparse(self.X):
            from scipy.sparse import vstack
            X = vstack(Xs)
        else:
            X = np.concatenate(Xs)
        obs = pd.concat([ad.obs for ad in adatas_to_concat])
        obsm = np.concatenate([ad.obsm for ad in adatas_to_concat])
        var = adatas_to_concat[0].var
        varm = adatas_to_concat[0].varm
        uns = adatas_to_concat[0].uns
        return AnnData(X, obs, var, uns, obsm, varm, filename=self.filename)

    concatenate.__doc__ = dedent("""\
        Concatenate along the observations axis after intersecting the variables names.

        The `.var`, `.varm`, and `.uns` attributes of the passed adatas are ignored.

        Parameters
        ----------
        adatas : :class:`~anndata.AnnData` or list of :class:`~anndata.AnnData`
            AnnData matrices to concatenate with.
        batch_key : `str` (default: 'batch')
            Add the batch annotation to `.obs` using this key.
        batch_categories : list (default: `range(len(adatas)+1)`)
            Use these as categories for the batch annotation.

        Returns
        -------
        adata : :class:`~anndata.AnnData`
            The concatenated AnnData, where `adata.obs['batch']` stores a
            categorical variable labeling the batch.

        Examples
        --------
        {example_concatenate}
        """).format(example_concatenate=_example_concatenate)

    def __contains__(self, key):
        raise AttributeError('AnnData has no attribute __contains__, don\'t check `in adata`.')

    def _check_dimensions(self, key=None):
        if key is None:
            key = {'obs', 'var', 'obsm', 'varm'}
        else:
            key = {key}
        if 'obs' in key and len(self._obs) != self._n_obs:
            raise ValueError('Observations annot. `obs` must have number of '
                             'rows of `X` ({}), but has {} rows.'
                             .format(self._n_obs, self._obs.shape[0]))
        if 'var' in key and len(self._var) != self._n_vars:
            raise ValueError('Variables annot. `var` must have number of '
                             'columns of `X` ({}), but has {} rows.'
                             .format(self._n_vars, self._var.shape[0]))
        if 'obsm' in key and len(self._obsm) != self._n_obs:
            raise ValueError('Observations annot. `obsm` must have number of '
                             'rows of `X` ({}), but has {} rows.'
                             .format(self._n_obs, self._obs.shape[0]))
        if 'varm' in key and len(self._varm) != self._n_vars:
            raise ValueError('Variables annot. `varm` must have number of '
                             'columns of `X` ({}), but has {} rows.'
                             .format(self._n_vars, self._var.shape[0]))

    def write(self, filename=None, compression='gzip',
              compression_opts=None):
        """Write `.h5ad`-formatted hdf5 file.

        Parameters
        ----------
        filename : `str`, optional (default: AnnData.filename)
            Filename of data file. Defaults to filename of object.
        compression : `None` or {'gzip', 'lzf'}, optional (default: `'gzip'`)
            See http://docs.h5py.org/en/latest/high/dataset.html.
        compression_opts : `int`, optional (default: `None`)
            See http://docs.h5py.org/en/latest/high/dataset.html.
        """
        if filename is None and self.filename is None:
            raise ValueError('Provide a filename!')
        if filename is None:
            filename = self.filename
        if self.filename is not None:
            X = self._read('X')
            adata = AnnData(X, self.obs, self.var, self.uns, self.obsm, self.varm)
        else:
            adata = self
        write_h5ad(filename, adata, compression, compression_opts)

    def write_csvs(self, dirname, skip_data=True):
        """Write annotation to `.csv` files.

        It is not possible to recover the full :class:`~anndata.AnnData` from the
        output of this function. Use :func:`~anndata.write` for this.

        Parameters
        ----------
        dirname : `str`
            Name of directory to which to export.
        skip_data : `bool`, optional (default: True)
             Skip the data matrix `.X`.
        """
        from .readwrite.write import write_csvs
        write_csvs(dirname, self, skip_data)

    def write_loom(self, filename):
        """Write `.loom`-formatted hdf5 file.

        Parameters
        ----------
        filename : `str`
            The filename.
        """
        from .readwrite.write import write_loom
        write_loom(filename, self)

    def _to_dict_fixed_width_arrays(self):
        """A dict of arrays that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        obs_rec, uns_obs = df_to_records_fixed_width(self._obs)
        var_rec, uns_var = df_to_records_fixed_width(self._var)
        d = {'obs': obs_rec,
             'var': var_rec,
             'obsm': self._obsm,
             'varm': self._varm}
        if self.filename is None: d['X'] = self._X
        else: d['X'] = self._read('X')
        return {**d, **self._uns, **uns_obs, **uns_var}

    def _from_dict(self, ddata):
        """Allows to construct an instance of AnnData from a dictionary.

        In particular, from a dict that has been written using
        ``AnnData._to_dict_fixed_width_arrays``.
        """
        ddata = ddata.copy()
        # data matrix
        if '_data' in ddata:
            X = ddata['_data']
            del ddata['_data']
        elif 'data' in ddata:
            X = ddata['data']
            del ddata['data']
        elif '_X' in ddata:
            X = ddata['_X']
            del ddata['_X']
        elif 'X' in ddata:
            X = ddata['X']
            del ddata['X']
        else:
            X = None
            
        # simple annotation
        if ('_obs' in ddata and isinstance(ddata['_obs'], (np.ndarray, pd.DataFrame))
            and '_var' in ddata and isinstance(ddata['_var'], (np.ndarray, pd.DataFrame))):
            obs = ddata['_obs']
            del ddata['_obs']
            var = ddata['_var']
            del ddata['_var']
        elif ('_smp' in ddata and isinstance(ddata['_smp'], (np.ndarray, pd.DataFrame))
              and '_var' in ddata and isinstance(ddata['_var'], (np.ndarray, pd.DataFrame))):
            obs = ddata['_smp']
            del ddata['_smp']
            var = ddata['_var']
            del ddata['_var']
        elif ('obs' in ddata and isinstance(ddata['obs'], (np.ndarray, pd.DataFrame))
              and 'var' in ddata and isinstance(ddata['var'], (np.ndarray, pd.DataFrame))):
            obs = ddata['obs']
            del ddata['obs']
            var = ddata['var']
            del ddata['var']
        elif ('smp' in ddata and isinstance(ddata['smp'], (np.ndarray, pd.DataFrame))
              and 'var' in ddata and isinstance(ddata['var'], (np.ndarray, pd.DataFrame))):
            obs = ddata['smp']
            del ddata['smp']
            var = ddata['var']
            del ddata['var']
        else:
            obs, var = OrderedDict(), OrderedDict()
            if 'row_names' in ddata:
                obs['obs_names'] = ddata['row_names']
                del ddata['row_names']
            elif 'obs_names' in ddata:
                obs['obs_names'] = ddata['obs_names']
                del ddata['obs_names']
            elif 'smp_names' in ddata:
                obs['obs_names'] = ddata['smp_names']
                del ddata['smp_names']
            if 'col_names' in ddata:
                var['var_names'] = ddata['col_names']
                del ddata['col_names']
            elif 'var_names' in ddata:
                var['var_names'] = ddata['var_names']
                del ddata['var_names']
            obs = {**obs, **ddata.get('row', {}), **ddata.get('obs', {})}
            var = {**var, **ddata.get('col', {}), **ddata.get('var', {})}
            for k in ['row', 'obs', 'col', 'var']:
                if k in ddata:
                    del ddata[k]

        # transform recarray to dataframe
        if isinstance(obs, np.ndarray) and isinstance(var, np.ndarray):
            from pandas.api.types import is_string_dtype
            from pandas import Index

            if 'obs_names' in obs.dtype.names:
                obs = pd.DataFrame.from_records(obs, index='obs_names')
            elif 'smp_names' in obs.dtype.names:
                obs = pd.DataFrame.from_records(obs, index='smp_names')
            elif 'row_names' in obs.dtype.names:
                obs = pd.DataFrame.from_records(obs, index='row_names')
            elif 'index' in obs.dtype.names:
                obs = pd.DataFrame.from_records(obs, index='index')
            else:
                obs = pd.DataFrame.from_records(obs)
            obs.index = obs.index.astype('U')
            # transform to unicode string
            # TODO: this is quite a hack
            for c in obs.columns:
                if is_string_dtype(obs[c]):
                    obs[c] = Index(obs[c]).astype('U').values

            if 'var_names' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='var_names')
            elif 'col_names' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='col_names')
            elif 'index' in var.dtype.names:
                var = pd.DataFrame.from_records(var, index='index')
            else:
                var = pd.DataFrame.from_records(var)
            var.index = var.index.astype('U')
            for c in var.columns:
                if is_string_dtype(var[c]):
                    var[c] = Index(var[c]).astype('U').values

            # these are the category fields
            k_to_delete = []
            for k in ddata.keys():
                if k.endswith('_categories'):
                    k_stripped = k.replace('_categories', '')
                    if k_stripped in obs:
                        obs[k_stripped] = pd.Categorical.from_codes(
                            codes=obs[k_stripped].values,
                            categories=ddata[k])
                    if k_stripped in var:
                        var[k_stripped] = pd.Categorical.from_codes(
                            codes=var[k_stripped].values,
                            categories=ddata[k])
                    k_to_delete.append(k)

            for k in k_to_delete:
                del ddata[k]

        # multicolumn annotation
        if 'obsm' in ddata:
            obsm = ddata['obsm']
            del ddata['obsm']
        elif '_obsm' in ddata:
            obsm = ddata['_obsm']
            del ddata['_obsm']
        elif 'smpm' in ddata:
            obsm = ddata['smpm']
            del ddata['smpm']
        elif '_smpm' in ddata:
            obsm = ddata['_smpm']
            del ddata['_smpm']
        else:
            obsm = None
        if 'varm' in ddata:
            varm = ddata['varm']
            del ddata['varm']
        elif '_varm' in ddata:
            varm = ddata['_varm']
            del ddata['_varm']
        else:
            varm = None

        # the remaining fields are the unstructured annotation
        uns = ddata

        return X, obs, var, uns, obsm, varm

    def _clean_up_old_format(self, uns):
        # multicolumn keys
        # all of the rest is only for backwards compat
        if uns and '_obs_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['_obs_keys_multicol'])
            del uns['_obs_keys_multicol']
        elif uns and 'obs_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['obs_keys_multicol'])
            del uns['obs_keys_multicol']
        elif uns and '_smp_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['_smp_keys_multicol'])
            del uns['_smp_keys_multicol']
        elif uns and 'smp_keys_multicol' in uns:
            _keys_multicol_obs = list(uns['smp_keys_multicol'])
            del uns['smp_keys_multicol']
        else:
            _keys_multicol_obs = []
        if uns and '_var_keys_multicol' in uns:
            _keys_multicol_var = list(uns['_var_keys_multicol'])
            del uns['_var_keys_multicol']
        elif uns and 'var_keys_multicol' in uns:
            _keys_multicol_var = list(uns['var_keys_multicol'])
            del uns['var_keys_multicol']
        else:
            _keys_multicol_var = []

        # now, for compat, fill the old multicolumn entries into obsm and varm
        # and remove them from obs and var
        for key in _keys_multicol_obs:
            self._obsm[key] = self._get_multicol_field_obs(key)
        for key in _keys_multicol_var:
            self._varm[key] = self._get_multicol_field_var(key)

    def _get_multicol_field_obs(self, key_multicol):
        return self._get_and_delete_multicol_field('obs', key_multicol)

    def _get_multicol_field_var(self, key_multicol):
        return self._get_and_delete_multicol_field('var', key_multicol)

    def _get_and_delete_multicol_field(self, a, key_multicol):
        keys = []
        for k in getattr(self, a).columns:
            if k.startswith(key_multicol):
                keys.append(k)
        values = getattr(self, a)[keys].values
        getattr(self, a).drop(keys, axis=1, inplace=True)
        return values


AnnData.__doc__ = dedent("""\
    An annotated data matrix.

    {main_narrative}

    Parameters
    ----------
    X : `np.ndarray`, `sp.spmatrix`, `np.ma.MaskedArray`
        A #observations × #variables data matrix. A view of the data is used if the
        data type matches, otherwise, a copy is made.
    obs : `pd.DataFrame`, `dict`, structured `np.ndarray` or `None`, optional (default: `None`)
        Key-indexed one-dimensional observation annotation of length #observations.
    var : `pd.DataFrame`, `dict`, structured `np.ndarray` or `None`, optional (default: `None`)
        Key-indexed one-dimensional variable annotation of length #variables.
    uns : `dict` or `None`, optional (default: `None`)
        Unstructured annotation for the whole dataset.
    obsm : structured `np.ndarray`, optional (default: `None`)
        Key-indexed multi-dimensional observation annotation of length #observations.
    varm : structured `np.ndarray`, optional (default: `None`)
        Key-indexed multi-dimensional observation annotation of length #observations.
    dtype : simple `np.dtype`, optional (default: 'float32')
        Data type used for storage.
    single_col : `bool`, optional (default: `False`)
        Interpret one-dimensional input array as column.

    See Also
    --------
    read
    read_csv
    read_excel
    read_hdf
    read_loom
    read_mtx
    read_text

    Notes
    -----
    Multi-dimensional annotations are stored in ``.obsm`` and ``.varm``.

    :class:`~anndata.AnnData` stores observations (samples) of variables
    (features) in the rows of a matrix. This is the convention of the modern
    classics of stats [Hastie09]_ and Machine Learning `(Murphy, 2012)
    <https://mitpress.mit.edu/books/machine-learning-0>`_, the convention of
    dataframes both in R and Python and the established stats and machine
    learning packages in Python (`statsmodels
    <http://www.statsmodels.org/stable/index.html>`_, `scikit-learn
    <http://scikit-learn.org/>`_). It is the opposite of the convention for
    storing genomic data.

    Examples
    --------
    A data matrix is flattened if either #observations (`n_obs`) or #variables
    (`n_vars`) is 1, so that Numpy's slicing behavior is reproduced::

        adata = AnnData(np.ones((2, 2)))
        adata[:, 0].X == adata.X[:, 0]

    AnnData objects can be concatenated via :func:`~anndata.AnnData.concatenate`.

    {example_concatenate}

    """).format(main_narrative=AnnData._main_narrative,
                example_concatenate=AnnData._example_concatenate)
