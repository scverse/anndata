"""Wraps `h5py <http://www.h5py.org/>`_ to handle sparse matrices.

:mod:`anndata.h5py` is based on and uses the conventions of
`h5sparse <https://github.com/appier/h5sparse>`_ by
`Appier Inc. <https://www.appier.com/>`_.
See the copyright and license note in the source code.

The design choices of :mod:`anndata.h5py`, however, are different. In
particular, :mod:`anndata.h5py` allows handling sparse and non-sparse data at
the same time. It achieves this by extending the functionality of :class:`File`
and :class:`Group` objects in the h5py API, and by providing a new
:class:`SparseDataset` object.

For examples and further information, see this `blog post <https://falexwolf.de/blog/171212_sparse_matrices_with_h5py>`_.

.. autosummary::
  :toctree: .

   File
   Group
   Dataset
   SparseDataset
"""
from .h5sparse import File, Group, SparseDataset, _load_h5_dataset_as_sparse
from h5py import Dataset, special_dtype

# Problem: the H5py intersphinx is broken, and only contains e.g. `Dataset` directly.
# So we can’t possibly link to it using :class:`Dataset`, since that will always find our version.
Dataset.__doc__ = """Equivalent to :class:`h5py.Dataset <h5py:Dataset>`."""
