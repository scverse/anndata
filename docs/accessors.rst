.. _accessors:

Accessors and paths
===================

.. module:: anndata.acc

:mod:`!anndata.acc` creates references to 1D and 2D arrays in an :class:`~anndata.AnnData` object.
You can use these to drive e.g. plotting or validation code.
For these purposes, they

#.  are easy to create:

    The central :attr:`A` object allows you to create
    :class:`AdRef` objects that reference arrays
    along one or two axes of an :class:`anndata.AnnData` object:

    >>> from anndata.acc import A
    >>> A[:, "gene-3"]
    A[:, 'gene-3']
    >>> type(A[:, "gene-3"])
    <class 'anndata.acc.AdRef'>

#.  introspectible

    :class:`AdRef`\ s have the :attr:`AdRef.axes`, :attr:`AdRef.idx`, and :attr:`AdRef.acc` attributes,
    allowing you to inspect all relevant properties.

    >>> A.var["symbol"].axes
    {'var'}
    >>> pc0 = A.obsm["pca"][:, 0]
    >>> pc0.idx
    0
    >>> pc0.acc
    A.obsm['pca']
    >>> pc0.acc.k
    'pca'

#.  convenient

    Want to reference multiple vectors from the same object?
    Pass a list of indices to the vector accessor:

    >>> A.obsp["connectivities"][:, ["cell0", "cell1"]]
    [A.obsp['connectivities'][:, 'cell0'], A.obsp['connectivities'][:, 'cell1']]

#.  extensible (see :`extending accessors`_)

API
---

Most importantly, there are

.. autodata:: A

..  autosummary::
    :toctree: generated/

    AdRef

The following classes are behind :attr:`AdRef.acc`,
and therefore useful in :ref:`matches <match>` or :func:`isinstance` checks:

..  autosummary::
    :toctree: generated/

    VecAcc
    MetaVecAcc
    LayerVecAcc
    MultiVecAcc
    GraphVecAcc

Finally, these classes are only useful for extending:

..  autosummary::
    :toctree: generated/

    AdAcc
    LayerAcc
    MultiAcc
    GraphAcc

.. _extending accessors:

Extending accessors
-------------------

.. attention:: TODO
