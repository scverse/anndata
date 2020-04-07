Concatenation
=============

.. note::

    This section is currently under construction. There will be more here in the future.

Merging `.uns`
--------------

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
