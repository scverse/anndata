# Write compatibility profiles

Every so often, `anndata` changes what it writes to disk:
it starts using an encoding that older versions can’t read,
or it starts rejecting output that violates the {doc}`on-disk spec <fileformat-prose>`.

Every write function takes a `compat` argument selecting which anndata version’s write behavior to target,
as a {class}`~anndata.WriteCompat` member or its string value ({data}`~anndata.WriteCompatStr`):

```python
adata.write_h5ad("old.h5ad", compat="0.12")
```

A profile is the oldest anndata version whose I/O behavior we promise to be compatible with.
Choosing an old one therefore does two things:

1. Encodings that version can’t read are disallowed.
2. Restrictions added after that version are relaxed.

## Profiles

```{list-table}
:header-rows: 1

* - Profile
  - Released
  - Effect
* - `"0.12"`
  - 2025-07-16
  - Relaxes both restrictions added in 0.13:
    3-dimensional {attr}`~anndata.AnnData.X`/{attr}`~anndata.AnnData.layers` ({issue}`2430`) and
    `/` in `h5ad` keys of `obs`, `var` and `uns` ({issue}`2039`) can be written again.
* - `"0.13"`
  - 2026-07-07
  - No relaxations. This is the default.
```

Both relaxations produce files that violate the spec – anndata wrote such files before 0.13
and still reads them, but other older implementations may not.
`/` in keys becomes a nested group on disk, which is why 0.13 banned it,
and `zarr` never allowed it in the first place.

To write files for anndata 0.11 or older, use anndata 0.13 or older:
0.12 is the oldest profile since e.g. reading a `zarr` store with those versions
requires `zarr` v2, which anndata ≥0.14 no longer writes.

## Checking before writing

{meth}`anndata.AnnData.unwriteable` takes the same `compat` argument,
so you can check an object against a profile without writing anything:

```python
if adata.unwriteable(compat="0.12"):
    ...
```

## The default profile

The default is {attr}`anndata.WriteCompat.DEFAULT`.

The default keeps lagging behind the newest profile by roughly a year,
so files stay readable by the anndata versions people are likely to have.
Pass the newest profile as `compat`
to opt into everything the installed anndata can write.
