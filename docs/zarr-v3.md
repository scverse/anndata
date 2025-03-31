# zarr-v3 Guide/Roadmap

`anndata` now uses the much improved {mod}`zarr` v3 package and also [allows writing of datasets in the v3 format](https://anndata.readthedocs.io/en/stable/generated/anndata.settings.html#anndata.settings.zarr_write_format), with the exception of structured arrays.  Users should notice a significant performance improvement, especially for cloud data, but also likely for local data as well.  Here is a quick guide on some of our learnings so far:

## Remote data

We now provide the {func}`anndata.experimental.read_lazy` feature for reading as much of the {class}`~anndata.AnnData` object as lazily as possible, using {mod}`dask` and {mod}`xarray`.  Please note that this feature is experimental and subject to change.  To enable this functionality in a performant and feature-complete way for remote data sources, we use [consolidated metadata](https://zarr.readthedocs.io/en/stable/user-guide/consolidated_metadata.html) on the `zarr` store (written by default).  Please note that this introduces consistency issues - if you update the structure of the underlying `zarr` store i.e., remove a column from `obs`, the consolidated metadata will no longer be valid.  Further, note that without consolidated metadata, we cannot guarantee your stored `AnnData` object will be fully readable.  And even if it is fully readable, it will almost certainly be much slower to read.

There are two ways of opening remote [`zarr` stores from the `zarr-python` package](https://zarr.readthedocs.io/en/stable/api/zarr/storage/index.html), `fsspec` and `obstore`, and both can be used with `read_lazy`.  [`obstore` claims to be more performant out-of-the-box](https://developmentseed.org/obstore/latest/performance), but notes that this claim has not been benchmarked with the `uvloop` event loop, which itself claims to be 2X more performant than the default event loop for `python`.

## Local data

Local data generally poses a different set of challenges.  First, write speeds can be somewhat slow and second, the creation of many small files on a file system can slow down a filesystem.  For the "many small files" problem, `zarr` has introduced [sharding in the v3 file format](https://zarr.readthedocs.io/en/stable/user-guide/performance.html#sharding). Sharding requires knowledge of the array element you are writing, though, and therefore you will need to use {func}`anndata.experimental.write_dispatched` to use sharding:

```python
import anndata as ad
from collections.abc import Mapping
from typing import Any
import zarr

ad.settings.zarr_write_format = 3 # Absolutely crucial! Sharding is only for the v3 file format!

def write_sharded(group: zarr.Group, adata: ad.AnnData):
    def callback(func: ad.experimental.Write, g: zarr.Group, k: str, elem: ad.typing.RWAble, dataset_kwargs: Mapping[str, Any], iospec: ad.experimental.IOSpec):
        if iospec.encoding_type in { "array" }:
            dataset_kwargs = { "shards": tuple(int(2 ** (16 / len(elem.shape))) for _ in elem.shape), **dataset_kwargs}
            dataset_kwargs["chunks"] = tuple(i // 2 for i in dataset_kwargs["shards"])
        elif iospec.encoding_type in { "csr_matrix", "csc_matrix" }:
            dataset_kwargs = { "shards": (2**16,), "chunks": (2**8, ), **dataset_kwargs }
        func(g, k, elem, dataset_kwargs=dataset_kwargs)
    return ad.experimental.write_dispatched(group, "/", adata, callback=callback)
```

However, `zarr-python` can be slow with sharding throughput as well as writing throughput.  Thus if you wish to [speed up](https://github.com/LDeakin/zarr_benchmarks) either writing, sharding, or both (or receive a modest speed-boost for reading), a [bridge to the `zarr` implementation in Rust](https://zarrs-python.readthedocs.io/en/latest/) can help with that:

```
uv pip install zarrs
```

```python
import zarr
import zarrs
zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})
```

However, this pipeline is not compatible with all types of zarr store, especially remote stores and there are limitations on where rust can give a performance boost for indexing.  We therefore recommend this pipeline for writing full datasets and reading contiguous regions of said written data.

## Codecs

The default `zarr-python` v3 codec for `v3 file-format` is no longer `blosc` but `zstd`.  While `zstd` is more widespread, you may find its performance to not meet your old expectations.  Therefore, we recommend passing in the [`BloscCodec`](https://zarr.readthedocs.io/en/stable/api/zarr/codecs/index.html#zarr.codecs.BloscCodec) if you wish to return to the old behavior.

There is currently a bug with `numcodecs` that prevents data written from other non-numcodecs `zstd` implementations from being read in by the default zarr pipeline (to which the above rust pipeline falls back if it cannot handle a datatype or indexing scheme, like `vlen-string`): https://github.com/zarr-developers/numcodecs/issues/424.  Thus is may be advisable to use `BloscCodec` with `zarr` v3 file format data if you wish to use the rust-accelerated pipeline until this issue is resolved.

The same issue with `zstd` applies to data that may eventually be written by the GPU `zstd` implementation (see below).

## GPU i/o

At the moment, it is unlikely your `anndata` i/o will work if you use [`zarr.enable_gpu`](https://zarr.readthedocs.io/en/stable/user-guide/gpu.html#reading-data-into-device-memory).  It's *possible* dense data i/o i.e., using {func}`anndata.io.read_elem` will work as expected, but this functionality is untested - sparse data, awkward arrays, and dataframes will not.  `kvikio` currently provides a [`GDS`-enable store](https://docs.rapids.ai/api/kvikio/nightly/api/#kvikio.zarr.GDSStore) although there are no working compressors at the moment exported from the `zarr-python` package (work is [underway for `Zstd`](https://github.com/zarr-developers/zarr-python/pull/2863)).

We anticipate enabling officially supporting this functionality officially for dense data, sparse data, and possibly awkward arrays in the next minor release, 0.13.

## `async`

At the moment, `anndata` exports no `async` functions.  However, `zarr-python` has a fully `async` API and provides its own event-loop so that users like `anndata` can interact with a synchronous API while still beenfitting from `zarr-python`'s asynchronous functionality under that API.  We anticipate providing `async` versions of {func}`anndata.io.read_elem` and {func}`anndata.experimental.read_dispatched` so that users can download data asynchronously without using the `zarr-python` event loop.  We also would like to create an asynchronous partial reader to enable iterative streaming of a dataset.
