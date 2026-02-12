.. _accessors:

Accessors and paths
===================

.. module:: anndata.acc

:mod:`!anndata.acc` provides :term:`accessor`\ s that create :term:`reference`\ s
to 1D and 2D arrays in :class:`~anndata.AnnData` objects.
You can use these to drive e.g. plotting or validation code.
For these purposes, they are

#.  easy to create:

    The central :attr:`A` object is an :term:`accessor` for the whole :class:`~anndata.AnnData` object,
    and allows you to create :class:`AdRef` objects,
    which are :term:`reference`\ s to arrays spanning one or two dimensions
    of an :class:`~anndata.AnnData` object (without being bound to a specific object):

    >>> from anndata.acc import A
    >>> A.X[:, "gene-3"]  # reference to `adata[:, "gene-3"].X` as 1D vector
    A.X[:, 'gene-3']
    >>> type(A.X[:, "gene-3"])
    <class 'anndata.acc.AdRef'>

    … and to use:

    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc3k_processed()

    E.g. to check if `adata.varm["PCs"]` has at least 30 columns:

    >>> A.varm["PCs"][:, 30] in adata
    True

    or to extract the referenced vector:

    >>> ref = A.obs["louvain"]
    >>> adata[ref].categories[:2]  # doctest: +ELLIPSIS
    Index(['CD4 T cells', 'CD14+ Monocytes'], dtype=...)

#.  introspectible:

    :class:`AdRef`\ s have the :attr:`AdRef.dims`, :attr:`AdRef.idx`, and :attr:`AdRef.acc` attributes,
    allowing you to inspect all relevant properties.

    >>> pc0 = A.obsm["pca"][:, 0]
    >>> pc0
    A.obsm['pca'][:, 0]
    >>> pc0.idx
    0
    >>> pc0.acc
    A.obsm['pca']
    >>> A.var["symbol"].dims
    {'var'}
    >>> pc0.acc.k
    'pca'

#.  convenient:

    Want to reference multiple vectors from the same object?
    Pass a list of indices to the vector accessor:

    >>> A.obsp["connectivities"][:, ["cell0", "cell1"]]
    [A.obsp['connectivities'][:, 'cell0'], A.obsp['connectivities'][:, 'cell1']]

#.  extensible: see `extending accessors`_.

API & Glossary
--------------

The central :term:`accessor` is :data:`!A`:

.. autodata:: A

See :class:`AdAcc` for examples of how to use it to create :term:`reference`\ s (:attr:`AdRef`\ s).

..  autosummary::
    :toctree: generated/
    :template: class-minimal

    AdAcc
    AdRef

..  glossary::

    reference
        An instance of :class:`AdRef`.
        References a 1D or 2D array in :class:`~anndata.AnnData` objects.
        It is independent of individual objects and can be inspected,
        checked for equality, used as mapping keys, or applied to concrete objects,
        e.g. via `ref in adata` or `adata[ref]`.

    accessor
        An instance of any of the `*Acc` classes, i.e. :class:`AdAcc`,
        or subclasses of :class:`MapAcc` or :class:`RefAcc`.
        Can be descended into via attribute access to get deeper accessors
        (e.g. `A` → `A.obs`) or references (e.g. `A.obs.index`, `A.obs["c"]`).

    reference accessor
        :class:`!RefAcc` subclasses directly create :term:`reference`\ s (:class:`AdRef` instances).
        They can be accessed from these references using the :attr:`AdRef.acc` attribute,
        and are therefore useful in :ref:`matches <match>` or :func:`isinstance` checks:

        ..  autosummary::
            :toctree: generated/
            :template: class-minimal

            RefAcc

        ..  list-table::
            :header-rows: 1

            -   - Class
                - Attributes
                - Examples
            -   - :class:`LayerAcc`
                - :attr:`LayerAcc.k`
                - `A.X[:, :]`, `A.layers["c"][:, "g0"]`
            -   - :class:`MetaAcc`
                - :attr:`MetaAcc.dim`
                - `A.obs["a"]`, `A.var["b"]`
            -   - :class:`MultiAcc`
                - :attr:`MultiAcc.dim`, :attr:`MultiAcc.k`
                - `A.obsm["d"][:, 2]`
            -   - :class:`GraphAcc`
                - :attr:`GraphAcc.dim`, :attr:`GraphAcc.k`
                - `A.obsp["e"][:, "c1"]`, `A.vbsp["e"]["g0", :]`

    mapping accessor
        :class:`MapAcc` subclasses can be indexed with a string to create :term:`reference accessor`\ s,
        e.g. `A.layers` or `A.obsm` are both :class:`MapAcc`\ s,
        while `A.layers["a"]` is a :class:`LayerAcc` and `A.obsm["b"]` is a :class:`MultiAcc`.
        :class:`!MapAcc`\ s are mostly useful for extending,
        but might be useful for APIs that need to refer to a :class:`~collections.abc.Mapping` of arrays:

        ..  autosummary::
            :toctree: generated/
            :template: class-minimal

            MapAcc
            LayerMapAcc
            MultiMapAcc
            GraphMapAcc

..  this indented code makes autosummary generate the pages
    ..  autosummary::
        :toctree: generated/
        :template: class-minimal

        LayerAcc
        MetaAcc
        MultiAcc
        GraphAcc

    ..  autosummary::
        :toctree: generated/

        Idx2D

..  toctree::
    :hidden:

    generated/anndata.acc.LayerAcc
    generated/anndata.acc.MetaAcc
    generated/anndata.acc.MultiAcc
    generated/anndata.acc.GraphAcc
    generated/anndata.acc.Idx2D

.. _extending accessors:

Extending accessors
-------------------

There are three layers of extensibility:

#.  subclassing :class:`RefAcc` and creating a new :class:`AdRef` instance for creating them:

    ..  code-block:: python

        from matplotlib import pyplot as plt
        from anndata.acc import AdAcc, AdRef

        class MplRef(AdRef, str):
            """Matplotlib will only treat strings as references, so we subclass `str`."""
            def __new__(cls, acc, idx) -> None:
                obj = str.__new__(cls, str(AdRef(acc, idx)))
                AdRef.__init__(obj, acc, idx)
                return obj

        A = AdAcc(ref_class=MplRef)

        adata = sc.datasets.pbmc3k_processed()
        plt.scatter(*A.obsm["X_umap"][:, [0, 1]], c=A.obs["n_counts"], data=adata)


#.  subclass one or more of the :term:`reference accessor`\ s, and create a new :class:`AdAcc` instance:

    >>> from anndata.acc import AdAcc, AdRef, MetaAcc
    >>>
    >>> class TwoDRef(AdRef):
    ...     """A reference able to refer to multiple metadata columns."""
    ...     ...
    >>>
    >>> class MyMetaAcc(MetaAcc):
    ...     def __getitem__(self, k):
    ...         if isinstance(k, list):
    ...             # override default behavior of returning a list of refs
    ...             return self.ref_class(self, k)
    ...         return super().__getitem__(k)
    >>>
    >>> A = AdAcc(ref_class=TwoDRef, meta_cls=MyMetaAcc)
    >>> A.obs[["a", "b"]]
    A.obs[['a', 'b']]

#.  subclass :class:`AdAcc` to add new accessors:

    >>> from dataclasses import dataclass, field
    >>> from anndata.acc import AdAcc, MetaAcc
    >>>
    >>> @dataclass(frozen=True)
    ... class EHRAcc(AdAcc):
    ...     tem: MetaAcc = field(init=False)
    ...     def __post_init__(self) -> None:
    ...         super().__post_init__()
    ...         tem = MetaAcc("tem", ref_class=self.ref_class)
    ...         object.__setattr__(self, "tem", tem)  # necessary because it’s frozen
    >>>
    >>> A = EHRAcc()
    >>> A.tem["visit_id"]
    A.tem['visit_id']
