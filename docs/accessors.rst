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
    >>> A[:, "gene-3"]  # reference to `adata[:, 3].X` as 1D vector
    A[:, 'gene-3']
    >>> type(A[:, "gene-3"])
    <class 'anndata.acc.AdRef'>

#.  introspectible

    :class:`AdRef`\ s have the :attr:`AdRef.axes`, :attr:`AdRef.idx`, and :attr:`AdRef.acc` attributes,
    allowing you to inspect all relevant properties.

    >>> pc0 = A.obsm["pca"][:, 0]
    >>> pc0  # reference to the 0th column of `adata.obsm['pca']`
    A.obsm['pca'][:, 0]
    >>> pc0.idx
    0
    >>> pc0.acc
    A.obsm['pca']
    >>> A.var["symbol"].axes
    {'var'}
    >>> pc0.acc.k
    'pca'

#.  convenient

    Want to reference multiple vectors from the same object?
    Pass a list of indices to the vector accessor:

    >>> A.obsp["connectivities"][:, ["cell0", "cell1"]]
    [A.obsp['connectivities'][:, 'cell0'], A.obsp['connectivities'][:, 'cell1']]

#.  extensible (see `extending accessors`_)

API
---

The central starting point is :data:`!A`:

.. autodata:: A

See :class:`!AdAcc` for examples of how to use it:

..  autosummary::
    :toctree: generated/
    :template: class-minimal

    AdAcc
    AdRef

The following classes are behind :attr:`AdRef.acc`,
and therefore useful in :ref:`matches <match>` or :func:`isinstance` checks:

.. _reference-accessors:

..  autosummary::
    :toctree: generated/
    :template: class-minimal

    RefAcc

- :class:`MetaAcc` with :attr:`MetaAcc.ax`
- :class:`LayerAcc` with :attr:`LayerAcc.k`
- :class:`MultiAcc` with :attr:`MultiAcc.ax` and :attr:`MultiAcc.k`
- :class:`GraphAcc` with :attr:`GraphAcc.ax` and :attr:`GraphAcc.k`

..  hidden
    ..  autosummary::
        :toctree: generated/
        :template: class-minimal

        MetaAcc
        LayerAcc
        MultiAcc
        GraphAcc

.. toctree::
   :hidden:

   generated/anndata.acc.MetaAcc
   generated/anndata.acc.LayerAcc
   generated/anndata.acc.MultiAcc
   generated/anndata.acc.GraphAcc

Finally, these classes are only useful for extending:

..  autosummary::
    :toctree: generated/
    :template: class-minimal

    LayerMapAcc
    MultiMapAcc
    GraphMapAcc

.. _extending accessors:

Extending accessors
-------------------

.. attention:: TODO
