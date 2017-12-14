"""Wraps `h5py <http://www.h5py.org/>`_ to handle sparse matrices.

:class:`anndata.h5py` is based on and uses the conventions of `h5sparse
<https://github.com/appier/h5sparse>`_ by `Appier
Inc. <https://www.appier.com/>`_. See the copyright and license note in the
source code.

The design choices of :class:`anndata.h5py`, however, are fundamentally
different. In particular, :class:`anndata.h5py` allows handling sparse and
non-sparse data at the same time. It achieves this by extending the
functionality of `File` and `Group` objects in the h5py API, and by providing a
new `SparseDataset` object.

For examples and further information, see this `blog post <https://falexwolf.de/blog/171212_sparse_matrices_with_h5py>`_.

.. autosummary::
  :toctree: .

   File
   Group
   Dataset
   SparseDataset
"""
from .h5sparse import File, Group, SparseDataset
from h5py import Dataset

Dataset.__doc__ = """Equivalent to `h5py.Dataset`."""
