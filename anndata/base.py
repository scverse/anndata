"""Main class and helper functions.
"""
import os
from enum import Enum
from collections import Mapping, Sequence
import numpy as np
from numpy import ma
import pandas as pd
from pandas.core.index import RangeIndex
from scipy import sparse
from scipy.sparse import issparse
from scipy.sparse.sputils import IndexMixin
from textwrap import dedent
from . import h5py


class StorageType(Enum):
    Array = np.ndarray
    Masked = ma.MaskedArray
    Sparse = sparse.spmatrix

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
            if issparse(X): X = X.toarray()
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
        # open in backed-mode
        return AnnData(filename=filename)
    else:
        # load everything into memory
        # TODO: do this with less of a hack...
        #       the "0" currently replaces the "self"
        d = AnnData._read_h5ad(0, filename)
        return AnnData(d)


class AnnDataFileManager():
    """Lightweight backing file manager for AnnData.
    """

    def __init__(self, adata, filename=None):
        self._adata = adata
        self._open(filename)

    def __repr__(self):
        if self._filename is None:
            return 'Backing file manager: no file is set.'
        else:
            return 'Backing file manager of file {}.'.format(self._filename)

    def _open(self, filename):
        """Init the file manager."""
        self._filename = filename
        if filename is None:
            self._file = None
        else:
            self._file = h5py.File(filename, 'r+')

    def close(self):
        """Close the backing file, remember filename, do *not* change to memory mode."""
        self._file.close()

    def _to_memory_mode(self):
        """Close the backing file, forget filename, *do* change to memory mode."""
        self.adata.__X = self.adata._read('X')[()]
        self._file.close()
        self._file = None
        self._filename = None

    def isopen(self):
        """State of backing file."""
        if self._file is None:
            return False
        try:
            # try accessing the id attribute to see if the file is open
            return bool(self.file._file.id)
        except:
            return False


class AnnData(IndexMixin):

    _BACKED_ATTRS = ['X']

    # backwards compat
    _H5_ALIASES = {
        'X': {'X', '_X', 'data', '_data'},
        'obs': {'obs', '_obs', 'smp', '_smp'},
        'var': {'var', '_var'},
        'obsm': {'obsm', '_obsm', 'smpm', '_smpm'},
        'varm': {'varm', '_varm'},
    }

    _H5_ALIASES_NAMES = {
        'obs': {'obs_names', 'smp_names', 'row_names', 'index'},
        'var': {'var_names', 'col_names', 'index'},
    }

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

    _doc = dedent("""\
        An annotated data matrix.

        {main_narrative}

        Parameters
        ----------
        X : `np.ndarray`, `sparse.spmatrix`, `np.ma.MaskedArray`
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

        """).format(main_narrative=_main_narrative,
                    example_concatenate=_example_concatenate)

    def __init__(self, X=None, obs=None, var=None, uns=None,
                 obsm=None, varm=None,
                 filename=None,
                 dtype='float32', single_col=False):

        # init from file
        if filename is None:
            self.file = AnnDataFileManager(self, None)
        else:
            if any((X, obs, var, uns, obsm, varm)):
                raise ValueError(
                    'If initializing from `filename`, '
                    'no further arguments may be passed.')
            self.file = AnnDataFileManager(self, filename)
            # will read from backing file
            # what is returned is, at this stage, a dict
            # that needs to be processed
            X = self._read_h5ad()

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
            if issparse(X) or isinstance(X, ma.MaskedArray):
                # TODO: maybe use view on data attribute of sparse matrix
                #       as in readwrite.read_10x_h5
                if X.dtype != np.dtype(dtype): X = X.astype(dtype)
            else:  # is np.ndarray
                X = X.astype(dtype, copy=False)
            # data matrix and shape
            self.__X, self._n_obs, self._n_vars = _fix_shapes(X)
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

    def __sizeof__(self):
        size = 0
        for attr in ['_X', '_obs', '_var', '_uns', '_obsm', '_varm']:
            s = getattr(self, attr).__sizeof__()
            size += s
        return size

    def _gen_repr(self, n_obs, n_vars):
        if self.isbacked():
            backed_at = 'backed at \'{}\''.format(self.filename)
        else:
            backed_at = ''
        descr = (
            'AnnData object with n_obs × n_vars = {} × {} {}'
            .format(n_obs, n_vars, backed_at))
        for attr in ['obs_keys', 'var_keys', 'uns_keys', 'obsm_keys', 'varm_keys']:
            keys = getattr(self, attr)()
            if len(keys) > 0:
                descr += '\n    {} = {}'.format(attr, keys)
        return descr

    def __repr__(self):
        return self._gen_repr(self.n_obs, self.n_vars)

    def _from_dict(self, ddata):
        """Allows to construct an instance of AnnData from a dictionary.

        Acts as interface for the communication with the hdf5 file.

        In particular, from a dict that has been written using
        ``AnnData._to_dict_fixed_width_arrays``.
        """
        # TODO: remove the copy at some point, check consequences first...
        ddata = ddata.copy()
        d_true_keys = {}

        for true_key, keys in AnnData._H5_ALIASES.items():
            for key in keys:
                if key in ddata:
                    d_true_keys[true_key] = ddata[key]
                    del ddata[key]
                    break
            else:
                d_true_keys[true_key] = None

        # transform recarray to dataframe
        from pandas.api.types import is_string_dtype
        from pandas import Index

        for true_key, keys in AnnData._H5_ALIASES_NAMES.items():
            if d_true_keys[true_key] is not None:
                for key in keys:
                    if key in d_true_keys[true_key].dtype.names:
                        d_true_keys[true_key] = pd.DataFrame.from_records(
                            d_true_keys[true_key], index=key)
                        break
                d_true_keys[true_key].index = d_true_keys[true_key].index.astype('U')
                # transform to unicode string
                # TODO: this is quite a hack
                for c in d_true_keys[true_key].columns:
                    if is_string_dtype(d_true_keys[true_key][c]):
                        d_true_keys[true_key][c] = Index(
                            d_true_keys[true_key][c]).astype('U').values

        # these are the category fields
        k_to_delete = []
        for k in ddata.keys():
            if k.endswith('_categories'):
                k_stripped = k.replace('_categories', '')
                if k_stripped in d_true_keys['obs']:
                    d_true_keys['obs'][k_stripped] = pd.Categorical.from_codes(
                        codes=d_true_keys['obs'][k_stripped].values,
                        categories=ddata[k])
                if k_stripped in d_true_keys['var']:
                    d_true_keys['var'][k_stripped] = pd.Categorical.from_codes(
                        codes=d_true_keys['var'][k_stripped].values,
                        categories=ddata[k])
                k_to_delete.append(k)

        for k in k_to_delete:
            del ddata[k]

        # assign the variables
        X = d_true_keys['X']
        obs = d_true_keys['obs']
        obsm = d_true_keys['obsm']
        var = d_true_keys['var']
        varm = d_true_keys['varm']
        # the remaining fields are the unstructured annotation
        uns = ddata

        return X, obs, var, uns, obsm, varm

    def _to_dict_fixed_width_arrays(self):
        """A dict of arrays that stores data and annotation.

        It is sufficient for reconstructing the object.
        """
        obs_rec, uns_obs = df_to_records_fixed_width(self._obs)
        var_rec, uns_var = df_to_records_fixed_width(self._var)
        d = {'X': self.__X,
             'obs': obs_rec,
             'var': var_rec,
             'obsm': self._obsm,
             'varm': self._varm}
        return {**d, **self._uns, **uns_obs, **uns_var}

    def _read_h5ad(self, filename=None):
        """Return a dict with arrays for initializing AnnData.

        Parameters
        ----------
        filename : `str` or `None`, optional (default: `None`)
            Defaults to the objects filename if `None`.
        """
        # we need to be able to call the function
        # without reference to self
        filename_was_none = False
        if filename is None:
            filename_was_none = True
            filename = self.filename

        def postprocess_reading(key, value):
            if value is None:
                return key, value
            if value.ndim == 1 and len(value) == 1:
                value = value[0]
            if value.dtype.kind == 'S':
                value = value.astype(str)
                # recover a dictionary that has been stored as a string
                if len(value) > 0:
                    if value[0] == '{' and value[-1] == '}': value = eval(value)
            # transform byte strings in recarrays to unicode strings
            # TODO: come up with a better way of solving this, see also below
            if (key not in AnnData._H5_ALIASES['obs']
                and key not in AnnData._H5_ALIASES['var']
                and not isinstance(value, dict) and value.dtype.names is not None):
                new_dtype = [((dt[0], 'U{}'.format(int(int(dt[1][2:])/4)))
                              if dt[1][1] == 'S' else dt) for dt in value.dtype.descr]
                value = value.astype(new_dtype)
            return key, value

        filename = str(filename)  # allow passing pathlib.Path objects
        d = {}
        if filename_was_none:
            f = self.file._file
        else:
            # open in editable mode to fix old file formats
            f = h5py.File(filename, 'r+')
        # TODO: this will still fail if trying to read a sparse matrix group
        for key in f.keys():
            if filename_was_none and key in AnnData._BACKED_ATTRS:
                # we cannot simply assign an h5py dataset here, as we need to
                # distinguish between dense and different sparse data types
                value = None
            else:
                # the '()' means 'load everything into memory' (by contrast, ':'
                # only works if not reading a scalar type)
                value = f[key][()]
            key, value = postprocess_reading(key, value)
            d[key] = value
        # backwards compat: save X with the correct name
        if 'X' not in d:
            for key in AnnData._H5_ALIASES['X']:
                if key in d:
                    del f[key]
                    f.create_dataset('X', data=d[key])
                    break
        # backwards compat: store sparse matrices properly
        csr_keys = [key.replace('_csr_data', '')
                    for key in d if '_csr_data' in key]
        for key in csr_keys:
            d = load_sparse_csr(d, key=key)
            # delete the old representation from the file
            del_sparse_matrix_keys(f, key)
            if key in AnnData._H5_ALIASES['X'] and key != 'X':
                d['X'] = d[key]
                del d[key]
                key = 'X'
            f.close()
            # store the new representation
            f = h5py.File(filename, 'r+')
            f.create_dataset(key, data=d[key])
            if filename_was_none is not None:
                d[key] = None
            f.close()
            f = h5py.File(filename, 'r+')
        if not filename_was_none:
            f.close()
        return d

    def _write_h5ad(self, filename, **kwargs):

        def preprocess_writing(value):
            if issparse(value):
                return value
            elif isinstance(value, dict):
                # hack for storing dicts
                value = np.array([str(value)])
            else:
                # make sure value is an array
                value = np.array(value)
                # hm, why that?
                if value.ndim == 0: value = np.array([value])
            # make sure string format is chosen correctly
            if value.dtype.kind == 'U': value = value.astype(np.string_)
            return value

        filename = str(filename)  # allow passing pathlib.Path objects
        if not filename.endswith(('.h5', '.h5ad')):
            raise ValueError('Filename needs to end with \'.h5ad\'.')
        if self.isbacked():
            # close so that we can reopen below
            self.file._file.close()
        # create directory if it doesn't exist
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        d = self._to_dict_fixed_width_arrays()
        # store dense values
        with h5py.File(filename, 'a') as f:
            for key, value in d.items():
                if (value is None
                    # ignore arrays with empty dtypes
                    or not value.dtype.descr):
                    continue
                value = preprocess_writing(value)
                try:
                    if key in set(f.keys()):
                        if (f[key].shape == value.shape
                            and f[key].dtype == value.dtype):
                            f[key][()] = value
                            continue
                        else:
                            del f[key]
                    f.create_dataset(key, data=value, **kwargs)
                except TypeError:
                    # try writing it as byte strings
                    try:
                        if value.dtype.names is None:
                            if key in set(f.keys()):
                                if (f[key].shape == value.shape
                                    and f[key].dtype == value.dtype):
                                    f[key][()] = value.astype('S')
                                    continue
                                else:
                                    del f[key]
                            f.create_dataset(key, data=value.astype('S'), **kwargs)
                        else:
                            new_dtype = [(dt[0], 'S{}'.format(int(dt[1][2:])*4))
                                         for dt in value.dtype.descr]
                            if key in set(f.keys()):
                                if (f[key].shape == value.shape
                                    and f[key].dtype == value.dtype):
                                    f[key][()] = value.astype(new_dtype)
                                    continue
                                else:
                                    del f[key]
                            f.create_dataset(
                                key, data=value.astype(new_dtype), **kwargs)
                    except Exception as e:
                        warnings.warn('Could not save field with key = "{}" '
                                      'to hdf5 file.'.format(key))
        if self.isbacked():
            self.file._file = h5py.File(filename, 'r+')

    def isbacked(self):
        """`True` if object is backed on disk, `False` otherwise."""
        return self.filename is not None

    @property
    def filename(self):
        """Change to backing mode by setting the filename of a `.h5ad` file.

        - Setting the filename writes the stored data to disk.
        - Setting the filename when the filename was previously another name
          moves the backing file from the previous file to the new file. If you
          want to copy the previous file, use copy(filename='new_filename').
        """
        return self.file._filename

    @filename.setter
    def filename(self, filename):
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
                    os.rename(self.filename, filename)
                else:
                    # do nothing
                    return
            else:
                # change from memory to backing-mode
                # write the content of self to disk
                self.write(filename)
            # open new file for accessing
            self.file._open(filename)
            # as the data is stored on disk, we can safely set self.__X to None
            self.__X = None

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
        return self.X

    # for backwards compat
    @data.setter
    def data(self, value):
        print('DEPRECATION WARNING: use attribute `.X` instead of `.data`, '
              '`.data` will be removed in the future.')
        self.X = value

    @property
    def X(self):
        """Data matrix of shape `n_obs` × `n_vars` (`np.ndarray`, `sparse.spmatrix`)."""
        return self._X

    @property
    def _X(self):
        # we need this additional layer for AnnDataView
        if self.isbacked(): return self._read_backed('X')
        else: return self.__X

    @X.setter
    def X(self, value):
        if not self.shape == value.shape:
            raise ValueError('Data matrix has wrong shape {}, need to be {}'
                             .format(value.shape, self.shape))
        self._X = value

    @_X.setter
    def _X(self, value):
        if self.isbacked(): self._write_backed('X', value)
        else: self.__X = value

    def _read_backed(self, attr):
        return self.file._file[attr]

    def _write_backed(self, attr, value):
        if (self.file._file[attr].shape == value.shape
            and not isinstance(self.file._file[attr], h5py.SparseDataset)
            and not issparse(value)):
            self.file._file[attr] = value
        else:
            del self.file._file[attr]
            self.file._file.create_dataset(attr, data=value)

    @property
    def n_smps(self):
        """Deprecated name for number of observations."""
        return self._n_obs

    @property
    def n_obs(self):
        """Number of observations."""
        return self._n_obs

    @property
    def n_vars(self):
        """Number of variables/features."""
        return self._n_vars

    @property
    def obs(self):
        """One-dimensional annotation of observations (`pd.DataFrame`)."""
        return self._obs

    @obs.setter
    def obs(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError('Can only assign pd.DataFrame.')
        if len(value) != self.n_obs:
            raise ValueError('Length does not match.')
        self._obs = value

    # for backwards compat
    @property
    def smp(self):
        """Deprecated access to one-dimensional annotation of observations (`pd.DataFrame`)."""
        return self.obs

    # for backwards compat
    @smp.setter
    def smp(self, value):
        self.obs = value

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
        """Names of observations (`.obs.index`)."""
        return self._obs.index

    # backwards compat
    @property
    def smp_names(self):
        """Deprecated access to names of observations."""
        return self._obs.index

    @obs_names.setter
    def obs_names(self, names):
        self._obs.index = names

    # backwards compat
    @smp_names.setter
    def smp_names(self, names):
        self._obs.index = names

    @property
    def var_names(self):
        """Names of variables (`.var.index`)."""
        return self._var.index

    @var_names.setter
    def var_names(self, names):
        self._var.index = names

    def obs_keys(self):
        """List keys of observation annotation `.obs`."""
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
                    raise IndexError(
                        'Name "{}" is not valid observation/variable name.'
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
            return slice(start, stop, step)
        elif isinstance(index, (int, str)):
            return name_idx(index)
        elif isinstance(index, (Sequence, np.ndarray)):
            return np.fromiter(map(name_idx, index), 'int64')
        else:
            raise IndexError('Unknown index {!r} of type {}'
                             .format(index, type(index)))

    # TODO: this is not quite complete...
    def __delitem__(self, index):
        obs, var = self._normalize_indices(index)
        # TODO: does this really work?
        if not self.isbacked():
            del self.__X[obs, var]
        else:
            X = self._read_backed('X')
            del X[obs, var]
            self._write_backed('X', X)
        if var == slice(None):
            del self._obs.iloc[obs, :]
        if obs == slice(None):
            del self._var.iloc[var, :]

    def __getitem__(self, index):
        """Returns a sliced view of the object."""
        return self._getitem_view(index)

    def _getitem_view(self, index):
        oidx, vidx = self._normalize_indices(index)
        return AnnDataView(self, oidx, vidx)

    def _getitem_copy(self, index):
        oidx, vidx = self._normalize_indices(index)
        if not self.isbacked(): X = self.__X[oidx, vidx]
        else: X = self._read_backed('X')[oidx, vidx]
        obs_new = self._obs.iloc[oidx]
        obsm_new = self._obsm[oidx]
        var_new = self._var.iloc[vidx]
        varm_new = self._varm[vidx]
        assert obs_new.shape[0] == X.shape[0], (oidx, obs_new)
        assert var_new.shape[0] == X.shape[1], (vidx, var_new)
        uns_new = self._uns.copy()
        self._slice_uns_sparse_matrices_inplace(uns_new)
        return AnnData(X, obs_new, var_new, uns_new, obsm_new, varm_new)

    def _slice_uns_sparse_matrices_inplace(self, uns, oidx):
        # slice sparse spatrices of n_obs × n_obs in self.uns
        if not (isinstance(oidx, slice) and
                oidx.start is None and oidx.step is None and oidx.stop is None):
            for k, v in uns.items():
                if isinstance(v, sparse.spmatrix) and v.shape == (self._n_obs, self._n_obs):
                    uns[k] = v.tocsc()[:, obs].tocsr()[oidx, :]

    def _inplace_subset_var(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[:, index], but inplace.
        """
        if not self.isbacked():
            self.__X = self.__X[:, index]
            self._n_vars = self.__X.shape[1]
        else:
            X = self._read_backed('X')
            X = X[:, index]
            self._n_vars = X.shape[1]
            self._write_backed('X', X)
        self._var = self._var.iloc[index]
        # TODO: the following should not be necessary!
        self._varm = BoundRecArr(self._varm[index], self, 'varm')
        return None

    def _inplace_subset_obs(self, index):
        """Inplace subsetting along variables dimension.

        Same as adata = adata[index, :], but inplace.
        """
        if not self.isbacked():
            self.__X = self.__X[index, :]
            self._n_obs = self.__X.shape[0]
        else:
            X = self._read_backed('X')
            X = X[index, :]
            self._n_obs = X.shape[0]
            self._write_backed('X', X)
        self._slice_uns_sparse_matrices_inplace(self._uns, index)
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
        if not self.isbacked():
            self.__X[obs, var] = val
        else:
            X = self._read_backed('X')
            X[obs, var] = val
            self._write_backed('X', X)

    def __len__(self):
        return self.shape[0]

    def transpose(self):
        """Transpose whole object.

        Data matrix is transposed, observations and variables are interchanged.
        """
        if not self.isbacked(): X = self.__X
        else: X = _read_backed('X')
        if sparse.isspmatrix_csr(X):
            return AnnData(X.T.tocsr(), self._var, self._obs, self._uns,
                           self._varm.flipped(), self._obsm.flipped(),
                           filename=self.filename)
        return AnnData(X.T, self._var, self._obs, self._uns,
                       self._varm.flipped(), self._obsm.flipped(),
                       filename=self.filename)

    T = property(transpose)

    def copy(self, filename=None):
        """Full copy, optionally on disk."""
        if filename is None:
            if not self.isbacked(): X = self.__X
            else: X = _read_backed('X')
            return AnnData(X.copy() if not self.isbacked() else X,
                           self._obs.copy(),
                           self._var.copy(), self._uns.copy(),
                           self._obsm.copy(), self._varm.copy())
        else:
            from shutil import copyfile
            copyfile(self.filename, filename)
            return AnnData(filename=filename)

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
        if issparse(self.X):
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

    def write(self, filename=None, compression='gzip', compression_opts=None):
        """Write `.h5ad`-formatted hdf5 file and close a potential backing file.

        Parameters
        ----------
        filename : `str`, optional (default: ``.filename``)
            Filename of data file. Defaults to backing file.
        compression : `None` or {'gzip', 'lzf'}, optional (default: `'gzip'`)
            See http://docs.h5py.org/en/latest/high/dataset.html.
        compression_opts : `int`, optional (default: `None`)
            See http://docs.h5py.org/en/latest/high/dataset.html.
        """
        if filename is None and not self.isbacked():
            raise ValueError('Provide a filename!')
        if filename is None:
            filename = self.filename
        self._write_h5ad(filename,
            compression=compression, compression_opts=compression_opts)
        if self.isbacked():
            self.file.close()

    def flush(self):
        """Write `.h5ad`-formatted hdf5 file and leave the backing file open."""
        self._write_h5ad(filename, compression='gzip', compression_opts=None)

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


AnnData.__doc__ = AnnData._doc


# let's not make this a subclass, but use composition:
# the view only has a subset of AnnData's capabilities
class AnnDataView():

    def __init__(self, adata, oidx, vidx):
        self._adata = adata
        self._oidx = oidx
        self._vidx = vidx
        if isinstance(oidx, slice):
            self.n_obs = len(self._adata.obs.index[oidx])
        elif isinstance(oidx, (int, np.int64)):
            self.n_obs = 1
        else:
            self.n_obs = len(oidx)
        if isinstance(vidx, slice):
            self.n_vars = len(self._adata.obs.index[vidx])
        elif isinstance(vidx, (int, np.int64)):
            self.n_vars = 1
        else:
            self.n_vars = len(vidx)

    def __repr__(self):
        return 'View of ' + self._adata._gen_repr(self.n_obs, self.n_vars)

    @property
    def X(self):
        return self._adata._X[self._oidx, self._vidx]

    @X.setter
    def X(self, value):
        self._adata._X[self._oidx, self._vidx] = value

    @property
    def obs(self):
        return self._adata._obs.iloc[self._oidx]

    @obs.setter
    def obs(self, value):
        self._adata._obs.iloc[self._oidx] = value

    @property
    def obsm(self):
        return self._adata._obsm[self._oidx]

    @obsm.setter
    def obsm(self, value):
        self._adata._obsm[self._oidx] = value

    @property
    def var(self):
        return self._adata._var.iloc[self._vidx]

    @var.setter
    def var(self, value):
        self._adata._var.iloc[self._vidx] = value

    @property
    def varm(self):
        return self._adata._varm[self._vidx]

    @varm.setter
    def varm(self, value):
        self._adata._varm[self._vidx] = value

    @property
    def uns(self):
        uns_new = self._uns.copy()
        self._adata._slice_uns_sparse_matrices_inplace(uns_new, self._oidx)
        return uns_new

    def copy(self):
        return self._adata._getitem_copy((self._oidx, self._vidx))

    def __getitem__(self, index):
        raise ValueError('Object is view and cannot be accessed with [].')

    def __setitem__(self, index, value):
        raise ValueError('Object is view and cannot be accessed with [].')


# all for backwards compat...
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
    del_sparse_matrix_keys(d)
    return d

def del_sparse_matrix_keys(mapping):
    del mapping[key_csr + '_data']
    del mapping[key_csr + '_indices']
    del mapping[key_csr + '_indptr']
    del mapping[key_csr + '_shape']
