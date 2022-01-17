.. role:: small
.. role:: smaller


On `master` :small:`the future`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: IO Specification

Internal handling of IO has been overhauled.
This should make it much easier to support new datatypes, use partial access, and use `AnnData` internally in other formats.

- Each element should be tagged with an `encoding_type` and `encoding_version`. See updated docs on the :doc:`file format <fileformat-prose>`
- Support for nullable integer and boolean data arrays. More data types to come!

.. rubric:: Features

- Compatibility with `h5ad` files written from Julia :pr:`569` :smaller:`I Kats`
- Many logging messages that should have been warnings are now warnings :pr:`650` :smaller:`I Virshup`
- Significantly more efficient :func:`anndata.read_umi_tools` :pr:`661` :smaller:`I Virshup`
- Fixed deepcopy of a copy of a view retaining sparse matrix view mixin type :pr:`670` :smaller:`M Klein`
- In many cases :attr:`~anndata.AnnData.X` can now be `None` :pr:`463` :smaller:`R Cannoodt` :pr:`677` :smaller:`I Virshup`. Remaining work is documented in :issue:`467`.
- Removed hard `xlrd` dependency :smaller:`I Virshup`

.. rubric:: Experimental Module



.. rubric:: Bug fixes

- Fixed issue where `.copy` was creating sparse matrices views when copying :pr:`670` :smaller:`michalk8`


0.7.8 :small:`9 November, 2021`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Bug fixes

- Re-include test helpers :pr:`641` :smaller:`I Virshup`


0.7.7 :small:`9 November, 2021`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Bug fixes

- Fixed propagation of import error when importing `write_zarr` but not all dependencies are installed :pr:`579` :smaller:`R Hillje`
- Fixed issue with `.uns` sub-dictionaries being referenced by copies :pr:`576` :smaller:`I Virshup`
- Fixed out-of-bounds integer indices not raising :class:`IndexError` :pr:`630` :smaller:`M Klein`
- Fixed backed `SparseDataset` indexing with scipy 1.7.2 :pr:`638` :smaller:`I Virshup`

.. rubric:: Development processes

- Use PEPs 621 (standardized project metadata), 631 (standardized dependencies), and 660 (standardized editable installs) :pr:`639` :smaller:`I Virshup`

0.7.6 :small:`11 April, 2021`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: New features

- Added :meth:`anndata.AnnData.to_memory` for returning an in memory object from a backed one :pr:`470` :pr:`542` :smaller:`V Bergen` :smaller:`I Virshup`
- :meth:`anndata.AnnData.write_loom` now writes `obs_names` and `var_names` using the `Index`'s `.name` attribute, if set :pr:`538` :smaller:`I Virshup`

.. rubric:: Bug fixes

- Fixed bug where `np.str_` column names errored at write time :pr:`457` :smaller:`I Virshup`
- Fixed "value.index does not match parent’s axis 0/1 names" error triggered when a data frame is stored in obsm/varm after obs_names/var_names is updated :pr:`461` :smaller:`G Eraslan`
- Fixed `adata.write_csvs` when `adata` is a view :pr:`462` :smaller:`I Virshup`
- Fixed null values being converted to strings when strings are converted to categorical :pr:`529` :smaller:`I Virshup`
- Fixed handling of compression key word arguments :pr:`536` :smaller:`I Virshup`
- Fixed copying a backed `AnnData` from changing which file the original object points at :pr:`533` :smaller:`ilia-kats`
- Fixed a bug where calling `AnnData.concatenate` an `AnnData` with no variables would error :pr:`537` :smaller:`I Virshup`

.. rubric:: Deprecations

- Passing positional arguments to :func:`anndata.read_loom` besides the path is now deprecated :pr:`538` :smaller:`I Virshup`
- :func:`anndata.read_loom` arguments `obsm_names` and `varm_names` are now deprecated in favour of `obsm_mapping` and `varm_mapping` :pr:`538` :smaller:`I Virshup`


0.7.5 :small:`12 November, 2020`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Functionality

- Added ipython tab completion and a useful return from `.keys` to `adata.uns` :pr:`415` :smaller:`I Virshup`

.. rubric:: Bug fixes

- Compatibility with `h5py>=3` strings :pr:`444` :smaller:`I Virshup`
- Allow `adata.raw = None`, as is documented :pr:`447` :smaller:`I Virshup`
- Fix warnings from pandas 1.1 :pr:`425` :smaller:`I Virshup`

0.7.4 :small:`10 July, 2020`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Concatenation overhaul :pr:`378` :smaller:`I Virshup`

- New function :func:`anndata.concat` for concatenating `AnnData` objects along either observations or variables
- New documentation section: :doc:`concatenation`

.. rubric:: Functionality

- AnnData object created from dataframes with sparse values will have sparse `.X` :pr:`395` :smaller:`I Virshup`

.. rubric:: Bug fixes

- Fixed error from `AnnData.concatenate` by bumping minimum versions of numpy and pandas :issue:`385`
- Fixed colors being incorrectly changed when `AnnData` object was subset :pr:`388`

0.7.3 :small:`20 May, 2020`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Bug fixes

- Fixed bug where graphs used too much memory when copying :pr:`381` :smaller:`I Virshup`

0.7.2 :small:`15 May, 2020`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
