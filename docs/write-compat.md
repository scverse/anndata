# Write compatibility profiles

Every so often, `anndata` starts using an on-disk encoding that older versions can’t read.

Every write function takes a `compat` argument selecting which anndata version’s write behavior to target,
as a {class}`~anndata.WriteCompat` member or its string value ({data}`~anndata.WriteCompatStr`):

```python
adata.write_h5ad("old.h5ad", compat="0.13")
```

A profile is the oldest anndata version whose I/O behavior we promise to be compatible with,
i.e. encodings that version can’t read are disallowed.

## Profiles

```{list-table}
:header-rows: 1

* - Profile
  - Released
  - Effect
* - `"0.13"`
  - 2026-07-07
  - This is the default.
* - `"0.14"`
  - unreleased
  - Enables the {ref}`sequence <sequence>` and {ref}`accessor <accessor>` encodings,
    i.e. writing heterogeneous, ragged and nested sequences in `uns` ({issue}`1979`)
    as well as {mod}`anndata.acc` accessors,
    and {ref}`escapes keys <escaped-keys>` that a file system may choke on, e.g. `/` and `:`.
    Older versions can’t read either.
```

## Checking before writing

{meth}`anndata.AnnData.unwriteable` takes the same `compat` argument,
so you can check an object against a profile without writing anything:

```python
if adata.unwriteable(compat="0.13"):
    ...
```

## The default profile

The default is {attr}`anndata.WriteCompat.DEFAULT`, later features are disabled.

The default keeps lagging behind the newest profile by roughly a year,
so files stay readable by the anndata versions people are likely to have.
Pass the newest profile as `compat`
to opt into everything the installed anndata can write.
