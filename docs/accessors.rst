.. _accessors:

Accessors and paths
===================

.. module:: anndata.acc

:mod:`!anndata.acc` creates references to 1D and 2D arrays in an :class:`~anndata.AnnData` object.
You can use these to drive e.g. plotting or validation code.
For these purposes, they are

#.  easy to create:

    The central :attr:`A` object is an accessor for the whole :class:`~anndata.AnnData` object,
    and allows you to create :class:`AdRef` objects that reference arrays
    along one or two dimensions of an :class:`~anndata.AnnData` object:

    >>> from anndata.acc import A
    >>> A[:, "gene-3"]  # reference to `adata[:, "gene-3"].X` as 1D vector
    A[:, 'gene-3']
    >>> type(A[:, "gene-3"])
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

API
---

The central starting point is :data:`!A`:

.. autodata:: A

See :class:`AdAcc` for examples of how to use it to create :attr:`AdRef`\ s.

..  autosummary::
    :toctree: generated/
    :template: class-minimal

    AdRef

.. _reference accessors:

Reference accessors
~~~~~~~~~~~~~~~~~~~

The following :class:`!RefAcc` subclasses correspond to the various key-value stores that anndata exposes.
They can be accessed from references using :attr:`AdRef.acc`,
and are therefore useful in :ref:`matches <match>` or :func:`isinstance` checks:

..  autosummary::
    :toctree: generated/
    :template: class-minimal

    RefAcc

..  list-table::
    :header-rows: 1

    - - Class
      - Attributes
      - Examples
    - - :class:`AdAcc`
      - is a :class:`LayerAcc` with `.k=None`
      - `A["c1", :]`, `A[:, "g1"]`, `A[:, :]`
    - - :class:`LayerAcc`
      - :attr:`LayerAcc.k`
      - `A.layers["c"][:, "g0"]`
    - - :class:`MetaAcc`
      - :attr:`MetaAcc.dim`
      - `A.obs["a"]`, `A.var["b"]`
    - - :class:`MultiAcc`
      - :attr:`MultiAcc.dim`, :attr:`MultiAcc.k`
      - `A.obsm["d"][:, 2]`
    - - :class:`GraphAcc`
      - :attr:`GraphAcc.dim`, :attr:`GraphAcc.k`
      - `A.obsp["e"][:, "c1"]`, `A.vbsp["e"]["g0", :]`

..  hidden
    ..  autosummary::
        :toctree: generated/
        :template: class-minimal

        AdAcc
        LayerAcc
        MetaAcc
        MultiAcc
        GraphAcc

    ..  autosummary::
        :toctree: generated/

        Idx2D

.. toctree::
   :hidden:

   generated/anndata.acc.AdAcc
   generated/anndata.acc.LayerAcc
   generated/anndata.acc.MetaAcc
   generated/anndata.acc.MultiAcc
   generated/anndata.acc.GraphAcc
   generated/anndata.acc.Idx2D

Mapping accessors
~~~~~~~~~~~~~~~~~

Finally, these classes are mostly useful for extending,
but might be useful for APIs that take a reference to a :class:`collections.abc.Mapping`
of arrays:

..  autosummary::
    :toctree: generated/
    :template: class-minimal

    MapAcc
    LayerMapAcc
    MultiMapAcc
    GraphMapAcc

.. _extending accessors:

Extending accessors
-------------------

There are three layers of extensibility:

#.  subclassing :class:`RefAcc` and creating a new :class:`AdRef` instance for creating them:

    ..  code-block:: python

        import my_plotting_library as pl

        class AdDim(RefAcc, pl.Dimension): ...
        A = AdAcc(ref_class=AdDim)

        pl.scatter(adata, A[:, "Actb"], color=A.obs["cell_type"])


#.  subclass one or more of the `reference accessors`_, and create a new :class:`AdAcc` instance:

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
    ...         tem = MetaAcc("tem", ref_class=self.ref_class)
    ...         object.__setattr__(self, "tem", tem)  # necessary because it’s frozen
    >>>
    >>> A = EHRAcc()
    >>> A.tem["visit_id"]
    A.tem['visit_id']
