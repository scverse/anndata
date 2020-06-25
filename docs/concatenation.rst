Concatenation
=============

.. note::

    This section is currently a draft.

With function :func:`~anndata.concat`, `AnnData` objects can be combined via a composition of two operations: concatenation and merging.

* Concatenation is when we keep all sub elements of each object, and stack these elements in an ordered way.
* Merging is combining a set of collections into one resulting collection which contains elements from the objects.

Concatenation
-------------

Let's start off with an example:

    >>> import scanpy as sc, anndata as ad, numpy as np, pandas as pd
    >>> from scipy import sparse
    >>> from anndata import AnnData
    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> pbmc
    AnnData object with n_obs × n_vars = 700 × 765
        obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'G2M_score', 'phase', 'louvain'
        var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
        uns: 'bulk_labels_colors', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
        obsm: 'X_pca', 'X_umap'
        varm: 'PCs'
        obsp: 'distances', 'connectivities'

If we split this object up by clusters of observations, then stack those subsets we'll obtain the same values – just ordered differently.

    >>> groups = pbmc.obs.groupby("louvain").indices
    >>> pbmc_concat = ad.concat([pbmc[inds] for inds in groups.values()], merge="same")
    >>> assert np.array_equal(pbmc.X, pbmc_concat[pbmc.obs_names].X)  # TODO: Make obs names not get renamed by default
    >>> pbmc_concat
    AnnData object with n_obs × n_vars = 700 × 765
        obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'G2M_score', 'phase', 'louvain'
        var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
        obsm: 'X_pca', 'X_umap'
        varm: 'PCs'
        obsp: 'distances', 'connectivities'

Note that we concatenated along the observations by default, and that all elements aligned to the observations were concatenated as well.

Inner and outer joins
~~~~~~~~~~~~~~~~~~~~~

When the variables present in the objects to be concatenated aren't exactly the same, you can choose to take either the intersection or union of these variables.
This is otherwise called taking the `"inner"` (intersection) or `"outer"` (union) join.
For example, given two anndata objects with differing variables:

    >>> a = AnnData(sparse.eye(3), var=pd.DataFrame(index=list("abc")))
    >>> b = AnnData(sparse.eye(2), var=pd.DataFrame(index=list("ba")))
    >>> ad.concat([a, b], join="inner").X.toarray()
    array([[1., 0.],
           [0., 1.],
           [0., 0.],
           [0., 1.],
           [1., 0.]], dtype=float32)
    >>> ad.concat([a, b], join="outer").X.toarray()
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.],
           [0., 1., 0.],
           [1., 0., 0.]], dtype=float32)

The join argument is used for any element which has both (1) an axis being concatenated and (2) has an axis not being concatenated.
When concatenating along the `obs` dimension, this means elements of `.X`, `obs`, `.layers`, and `.obsm` will be affected by the choice of `join`.

To demonstrate this, let's say we're trying to combine a droplet based experiment with a spatial one.
When building a joint anndata object, we would still like to store the coordinates for the spatial samples.

    >>> coords = np.hstack([np.repeat(np.arange(10), 10), np.tile(np.arange(10), 10)]).T
    >>> spatial = AnnData(
    ...     sparse.random(5000, 10000, format="csr"), 
    ...     obsm={"coords": np.random.randn(5000, 2)}
    ... )
    >>> droplet = AnnData(sparse.random(5000, 10000, format="csr"))
    >>> combined = ad.concat([spatial, droplet], join="outer")
    >>> sc.pl.embedding(combined, "coords")  # doctest: +SKIP

.. TODO: Get the above plot to show up

Merging
-------

We also provide control over how elements which aren't aligned to the axis being concatenated are combined.
We'll refer to these other axis as the alternative axes here.
Selecting the elements aligned to alternative axes is done through merging, which we provide a few strategies for:

* The default, where no elements aligned to alternative axes are present in the result object
* `"same"`: Elements which are the same in each of the objects.
* `"unique"`: Elements for which there is only one possible value.
* `"first"`: The first element seen at each from each position.
* `"only"`: Elements which show up in only one of the objects.

These are similar to `xarray`'s compat_ argument.

.. _compat: http://xarray.pydata.org/en/stable/generated/xarray.merge.html#xarray.merge

We'll show how this works with a few examples for alternative axes, and then discuss how they work for `.uns` in more detail.

    >>> import scanpy as sc
    >>> blobs = sc.datasets.blobs(n_variables=30, n_centers=5)
    >>> sc.pp.pca(blobs)
    >>> blobs
    AnnData object with n_obs × n_vars = 640 × 30
        obs: 'blobs'
        uns: 'pca'
        obsm: 'X_pca'
        varm: 'PCs'

Now we will split this object by the present categories and recombine it to illustrate different merge strategies.

    >>> adatas = []
    >>> for group, idx in blobs.obs.groupby("blobs").indices.items():
    ...     sub_adata = blobs[idx].copy()
    ...     sub_adata.obsm["qc"], sub_adata.varm[f"{group}_qc"] = sc.pp.calculate_qc_metrics(
    ...         sub_adata, percent_top=(), inplace=False, log1p=False
    ...     )
    ...     adatas.append(sub_adata)
    >>> adatas[0]
    AnnData object with n_obs × n_vars = 128 × 30
        obs: 'blobs'
        uns: 'pca'
        obsm: 'X_pca', 'qc'
        varm: 'PCs', '0_qc'

`adatas` is now a list of datasets with disjoint sets of observations and a common set of variables.
Each object has had QC metrics computed, with observation-wise metrics stored under `"qc"` in `.obsm`, and variable-wise metrics stored with a unique key for each subset.
Taking a look at how this effects concatenation:

    >>> ad.concat(adatas)
    AnnData object with n_obs × n_vars = 640 × 30
        obs: 'blobs'
        obsm: 'X_pca', 'qc'
    >>> ad.concat(adatas, merge="same")
    AnnData object with n_obs × n_vars = 640 × 30
        obs: 'blobs'
        obsm: 'X_pca', 'qc'
        varm: 'PCs'
    >>> ad.concat(adatas, merge="unique")
    AnnData object with n_obs × n_vars = 640 × 30
        obs: 'blobs'
        obsm: 'X_pca', 'qc'
        varm: 'PCs', '0_qc', '1_qc', '2_qc', '3_qc', '4_qc'

Note that comparisons are made after indices are aligned.
That is, if the objects only share a subset of indices on the alternative axis, it's only required that values for those indices match when using a strategy like `"same"`.

    >>> a = AnnData(
    ...     sparse.eye(3),
    ...     var=pd.DataFrame({"nums": [1, 2, 3]}, index=list("abc"))
    ... )
    >>> b = AnnData(
    ...     sparse.eye(2),
    ...     var=pd.DataFrame({"nums": [2, 1]}, index=list("ba"))
    ... )
    >>> ad.concat([a, b], merge="same").var
       nums
    a     1
    b     2


Merging `.uns`
~~~~~~~~~~~~~~

We use the same set of strategies for merging `uns` as we do for entries aligned to an axis, but these strategies are applied recursively.
This is a little abstract, so we'll look at some examples of this. Here's our setup:

    >>> from anndata import AnnData
    >>> import numpy as np
    >>> a = AnnData(np.zeros((10, 10)), uns={"a": 1, "b": 2, "c": {"c.a": 3, "c.b": 4}})
    >>> b = AnnData(np.zeros((10, 10)), uns={"a": 1, "b": 3, "c": {"c.b": 4}})
    >>> c = AnnData(np.zeros((10, 10)), uns={"a": 1, "b": 4, "c": {"c.a": 3, "c.b": 4, "c.c": 5}})

For quick reference, these are the results from each of the merge strategies.
These are discussed in more depth below:

===========  =======================================================
`uns_merge`  Result
===========  =======================================================
`None`       `{}`
`"same"`     `{"a": 1, "c": {"c.b": 4}}`
`"unique"`   `{"a": 1, "c": {"c.b": 4, "c.c": 5, "c.a": 3}}`
`"only"`     `{"c": {"c.c": 5}}`
`"first"`    `{"a": 1, "b": 2, "c": {"c.b": 4, "c.c": 5, "c.a": 3}}`
===========  =======================================================

The default returns a fairly obvious result:

    >>> ad.concat([a, b, c]).uns == {}
    True

But let's take a look at the others in a bit more depth. Here, we'll be wrapping the output data in a `dict` for simplicity of the return value.

    >>> dict(ad.concat([a, b, c], uns_merge="same").uns)
    {'a': 1, 'c': {'c.b': 4}}

Here only the values for `uns["a"]` and `uns["c"]["c.b"]` were exactly the same, so only they were kept.
`uns["b"]` has a number of values and neither `uns["c"]["c.a"]` or `uns["c"]["c.b"]` appears in each `uns`.

A key feature to note is that comparisons are aware of the nested structure of `uns` and will be applied at any depth.
This is why `uns["c"]["c.b"]` was kept.

Merging `uns` in this way can be useful when there is some shared data between the objects being concatenated.
For example, if each was put through the same pipeline with the same parameters, those parameters used would still be present in the resulting object.

Now let's look at the behaviour of `unique`:

    >>> dict(ad.concat([a, b, c], uns_merge="unique").uns)
    {'a': 1, 'c': {'c.a': 3, 'c.b': 4, 'c.c': 5}}

The results here are a super-set of those from `"same"`. Note that there was only one possible value at each position in the resulting mapping.
That is, there were not alternative values present for `uns["c"]["c.c"]` even though it appeared only once.

This can be useful when the object's were both run through the same pipeline but contain specific metadata per object.
An example of this would be a spatial dataset, where the images are stored in `uns`.

    >>> dict(ad.concat([a, b, c], uns_merge="only").uns)
    {'c': {'c.c': 5}}

`uns["c"]["c.c"]` is the only value that is kept, since it is the only one which was specified in only one `uns`.

    >>> dict(ad.concat([a, b, c], uns_merge="first").uns)
    {'a': 1, 'b': 2, 'c': {'c.a': 3, 'c.b': 4, 'c.c': 5}}
 
In this case, the result has the union of the keys from all the starting dictionaries.
The value is taken from the first object to have a value at this key.
