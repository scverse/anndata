.. _accessors:

Accessors and paths
===================

.. module:: anndata.acc

:mod:`!anndata.acc` creates references to 1D and 2D arrays in an :class:`~anndata.AnnData` object.
You can use these to drive e.g. plotting or validation code.
For these purposes, they

#.  are easy to create:

    The central :attr:`A` object allows you to create
    :class:`AdPath` objects that reference arrays
    along one or two axes of an :class:`anndata.AnnData` object:

    >>> from anndata.acc import A
    >>> A[:, "gene-3"]
    A[:, 'gene-3']
    >>> type(A[:, "gene-3"])
    <class 'anndata.acc.AdPath'>

#.  introspectible

    :class:`AdPath`\ s have the :attr:`AdPath.axes`, :attr:`AdPath.idx`, and :attr:`AdPath.acc` attributes,
    allowing you to inspect all relevant properties.

    >>> A.var["symbol"].axes
    {'var'}
    >>> pc0 = A.obsm["pca"][:, 0]
    >>> pc0.idx
    (slice(None, None, None), 0)
    >>> pc0.acc
    A.obsm['pca']
    >>> pc0.acc.k
    'pca'

#.  convenient

    Want to reference multiple vectors from the same object?
    Pass a list of indices to the vector accessor:

    >>> A.obsp["connectivities"][:, ["cell0", "cell1"]]
    [A.obsp['connectivities'][:, 'cell0'], A.obsp['connectivities'][:, 'cell1']]

#. extensible (see :`extending accessors`_)

API
---

Most importantly, there are

.. autodata:: A

..  autosummary::
    :toctree: generated/

    anndata.acc.AdPath

The following classes are behind :attr:`AdPath.acc`,
and therefore useful in :ref:`matches <match>` or :func:`isinstance` checks:

..  autosummary::
    :toctree: generated/

    anndata.acc.VecAcc
    anndata.acc.LayerVecAcc
    anndata.acc.MultiVecAcc
    anndata.acc.GraphVecAcc

Finally, these classes are only useful for extending:

..  autosummary::
    :toctree: generated/

    anndata.acc.AdAcc
    anndata.acc.MultiAcc
    anndata.acc.GraphAcc

.. _extending accessors:

Extending accessors
-------------------

.. attention:: TODO
