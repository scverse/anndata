Concatenation
=============

.. note::

    This section is currently under construction. There will be more here in the future.

`AnnData` objects can be combined via a composition of two operations: concatenation and merging.
Concatenation is when we will keep all sub elements of each object to be combined, and stack these elements in an ordered way.
Merging is combining a set of collections into one resulting collection which contains elements from the inputs.

AnnData objects are combined using the function `ad.concat`. Here we'll look at how this function combines anndata objects using concatenation and merging.

Concatenation
-------------

This is essentially stacking anndata objects. Let's start off with an example:

    >>> import scanpy as sc, anndata as ad, numpy as np, pandas as pd
    >>> from anndata._core.merge import concat
    >>> from anndata import AnnData
    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> pbmc
    AnnData object with n_obs × n_vars = 700 × 765
        obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'G2M_score', 'phase', 'louvain'
        var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
        uns: 'louvain', 'neighbors', 'pca', 'rank_genes_groups'
        obsm: 'X_pca', 'X_umap'
        varm: 'PCs'
        obsp: 'distances', 'connectivities'

We can split this object up by the clusters of it's observations, then stack it again and we should the same values, just ordered differently.

    >>> groups = pbmc.obs.groupby("louvain").indices
    >>> pbmc_concat = ad.concat([pbmc[inds] for inds in groups.values()], index_unique=None)
    >>> assert np.array_equal(pbmc.X, pbmc_concat[pbmc.obs_names].X)  # TODO: Make obs names not get renamed by default
    >>> pbmc_concat
    AnnData object with n_obs × n_vars = 700 × 765
        obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'G2M_score', 'phase', 'louvain', 'batch'
        var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
        obsm: 'X_pca', 'X_umap'

Note that we concatenated along the observations by default, and that all elements aligned to the observations were concatenated as well.

Inner and outer joins
~~~~~~~~~~~~~~~~~~~~~

When the variables present in the objects to be concatnated aren't exactly the same, you can choose to take either the intersection or union of these variable.
This is otherwise called taking the `"inner"` (intersection) or `"outer"` (union) join.
For example, given two anndata objects with differing variables:

    >>> a = AnnData(sparse.eye(3), var=pd.DataFrame(index=list("abc")))
    >>> b = AnnData(sparse.eye(2), var=pd.DataFrame(index=list("ba")))
    >>> concat([a, b], join="inner").X.toarray()
    array([[1., 0.],
           [0., 1.],
           [0., 0.],
           [0., 1.],
           [1., 0.]], dtype=float32)
    >>> concat([a, b], join="outer").X.toarray()
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
            sparse.random(5000, 10000, format="csr"), 
            obsm={"coords": np.random.randn(5000, 2)}
        )
    >>> droplet = AnnData(sparse.random(5000, 10000, format="csr"))
    >>> combined = concat([spatial, droplet], join="outer")
    >>> sc.pl.embedding(combined, "coords")  # How to get this plot into the docs?

Merging
-------

Merging `.uns`
~~~~~~~~~~~~~~

When concatenating `AnnData` objects there are a number of options available for how the entries of `uns` should be merged.
These are:

* The default, creates a new empty `dict` for `uns`.
* `"same"`: Elements which are the same in each of the objects.
* `"unique"`: Elements for which there is only one possible value.
* `"first"`: The first element seen at each from each position.
* `"only"`: Elements which show up in only one of the objects.

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

    >>> a.concatenate([b, c]).uns == {}
    True

But let's take a look at the others in a bit more depth. Here, we'll be wrapping the output data in a `dict` for simplicity of the return value.

    >>> dict(a.concatenate([b, c], uns_merge="same").uns)
    {"a": 1, "c": {"c.b": 4}}

Here only the values for `uns["a"]` and `uns["c"]["c.b"]` were exactly the same, so only they were kept.
`uns["b"]` has a number of values and neither `uns["c"]["c.a"]` or `uns["c"]["c.b"]` appears in each `uns`.

A key feature to note is that comparisons are aware of the nested structure of `uns` and will be applied at any depth.
This is why `uns["c"]["c.b"]` was kept.

Merging `uns` in this way can be useful when there is some shared data between the objects being concatenated.
For example, if each was put through the same pipeline with the same parameters, those parameters used would still be present in the resulting object.

Now let's look at the behaviour of `unique`:

    >>> dict(a.concatenate([b, c], uns_merge="unique").uns)
    {"a": 1, "c": {"c.a": 3, "c.b": 4, "c.c": 5}}

The results here are a super-set of those from `"same"`. Note that there was only one possible value at each position in the resulting mapping.
That is, there were not alternative values present for `uns["c"]["c.c"]` even though it appeared only once.

This can be useful when the object's were both run through the same pipeline but contain specific metadata per object.
An example of this would be a spatial dataset, where the images are stored in `uns`.

    >>> dict(a.concatenate([b, c], uns_merge="only").uns)
    {"c": {"c.c": 5}}

`uns["c"]["c.c"]` is the only value that is kept, since it is the only one which was specified in only one `uns`.

    >>> dict(a.concatenate([b, c], uns_merge="first").uns)
    {"a": 1, "b": 2, "c": {"c.b": 4, "c.c": 5, "c.a": 3}}
 
In this case, the result has the union of the keys from all the starting dictionaries.
The value is taken from the first object to have a value at this key.
