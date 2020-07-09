.. role:: small
.. role:: smaller

0.7.4 :small:`2020-07-10`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Concatenation overhaul :pr:`378` :smaller:`I Virshup`

- New function :func:`anndata.concat` used for concatenating `AnnData` objects along either observations or variables
- New documentation section: :doc:`concatenation`

.. rubric:: Functionality

- AnnData object created from dataframes with sparse values will have sparse `.X` :pr:`395` :smaller:`I Virshup`

.. rubric:: Bug fixes

- Fixed error from `AnnData.concatenate` by bumping minimum versions of numpy and pandas :issue:`385`
- Fixed colors being incorrectly changed when `AnnData` object was subset :pr:`388`

0.7.3 :small:`2020-05-20`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Bug fixes

- Fixed bug where graphs used too much memory when copying :pr:`381` :smaller:`I Virshup`

0.7.2 :small:`2020-05-15`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Concatenation overhaul :smaller:`I Virshup`

- Elements of `uns` can now be merged, see :pr:`350`
- Outer joins now work for `layers` and `obsm`, see :pr:`352`
- Fill value for outer joins can now be specified
- Expect improvments in performance, see :issue:`303`

.. rubric:: Functionality

- :attr:`~anndata.AnnData.obsp` and :attr:`~anndata.AnnData.varp` can now be transposed :pr:`370` :smaller:`A Wolf`
- :meth:`~anndata.AnnData.obs_names_make_unique` is now better at making values unique, and will warn if ambiguities arise :pr:`345` :smaller:`M Weiden`
- :attr:`~anndata.AnnData.obsp` is now preferred for storing pairwise relationships between observations. In practice, this means there will be deprecation warnings and reformatting applied to objects which stored connectivities under `uns["neighbors"]`. Square matrices in :attr:`~anndata.AnnData.uns` will no longer be sliced (use `.{obs,var}p` instead). :pr:`337` :smaller:`I Virshup`
- :class:`~anndata.ImplicitModificationWarning` is now exported :pr:`315` :smaller:`P Angerer`
- Better support for :class:`~numpy.ndarray` subclasses stored in `AnnData` objects :pr:`335` :smaller:`michalk8`

.. rubric:: Bug fixes

- Fixed inplace modification of :class:`~pandas.Index` objects by the make unique function :pr:`348` :smaller:`I Virshup`
- Passing ambiguous keys to :meth:`~anndata.AnnData.obs_vector` and :meth:`~anndata.AnnData.var_vector` now throws errors :pr:`340` :smaller:`I Virshup`
- Fix instantiating :class:`~anndata.AnnData` objects from :class:`~pandas.DataFrame` :pr:`316` :smaller:`P Angerer`
- Fixed indexing into `AnnData` objects with arrays like `adata[adata[:, gene].X > 0]` :pr:`332` :smaller:`I Virshup`
- Fixed type of version :pr:`315` :smaller:`P Angerer`
- Fixed deprecated import from :mod:`pandas` :pr:`319` :smaller:`P Angerer`
