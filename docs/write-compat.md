# Write compatibility profiles

Every so often, the python package `anndata` starts using an on-disk encoding that older readers won’t be able to read,
whether that’s an older python `anndata` or an implementation in another language.

Every write function takes a `compat` argument selecting which anndata version’s write behavior to target,
as a {class}`~anndata.WriteCompat` member or its string value ({data}`~anndata.types.WriteCompatStr`):

```python
adata.write_h5ad("old.h5ad", compat="0.13")
```

A profile is the oldest anndata version whose I/O behavior we promise to be compatible with,
i.e. encodings that version can’t read are disallowed.
New encodings only land in minor releases, so e.g. `"0.13"` files are readable by every 0.13.x release.

This is different from an element’s `encoding-version` attribute (see {doc}`fileformat-prose`):
that versions the on-disk format of a single element type,
while a profile is a promise about the whole file: which `anndata` version can read it.

Profiles are named after versions of the python package,
but they describe sets of encodings, so they apply to other implementations too:
a reader that supports every encoding that python `anndata` X.Y can read
can read all files written with profile `"X.Y"`.

## Profiles

```{list-table}
:header-rows: 1

* - Profile
  - Released
  - Effect
* - `"0.13"`
  - 2026-07-07
  - This is the default.
```

## Checking before writing

{meth}`anndata.AnnData.unwriteable` takes the same `compat` argument,
so you can check an object against a profile without writing anything:

```python
if adata.unwriteable(compat="0.13"):
    ...
```

## The default profile

The default is {attr}`anndata.WriteCompat.DEFAULT`.

The default keeps lagging behind the newest profile by roughly a year,
so files stay readable by the anndata versions people are likely to have.
Pass the newest profile as `compat`
to opt into everything the installed anndata can write.
