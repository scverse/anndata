"""Wraps `h5py <http://www.h5py.org/>`_ to handle sparse matrices.

:mod:`anndata.mapping` is based on and uses the conventions of
`h5sparse <https://github.com/appier/h5sparse>`_ by
`Appier Inc. <https://www.appier.com/>`_.
See the copyright and license note in the source code.

The design choices of :mod:`anndata.mapping`, however, are different. In
particular, :mod:`anndata.h5py` allows handling sparse and non-sparse data at
the same time. It achieves this by extending the functionality of :class:`File`
and :class:`Group` objects in the h5py API, and by providing a new
:class:`SparseDataset` object.

For examples and further information, see this `blog post <https://falexwolf.de/blog/171212_sparse_matrices_with_h5py>`_.

.. autosummary::
  :toctree: .

   Group
   Dataset
   SparseDataset
"""
from .sparse_mapping import get_memory_class, add_mapping_impl
from .sparse_mapping import Group, SparseDataset

from .mapping_reader import MappingReader, report_key_on_error, AnnDataReadError
