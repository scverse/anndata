# zarr-v3 Guide/Roadmap

`anndata` now uses the much improved {mod}`zarr` v3 package and also allows writing of datasets in the v3 format via {attr}`anndata.settings.zarr_write_format`, with the exception of structured arrays.
Users should notice a significant performance improvement, especially for cloud data, but also likely for local data as well.
Here is a quick guide on some of our learnings so far:

## Remote data

We now provide the {func}`anndata.experimental.read_lazy` feature for reading as much of the {class}`~anndata.AnnData` object as lazily as possible, using `dask` and {mod}`xarray`.
Please note that this feature is experimental and subject to change.
To enable this functionality in a performant and feature-complete way for remote data sources, we use {doc}`zarr:user-guide/consolidated_metadata` on the `zarr` store (written by default).
Please note that this introduces consistency issues – if you update the structure of the underlying `zarr` store i.e., remove a column from `obs`, the consolidated metadata will no longer be valid.
Further, note that without consolidated metadata, we cannot guarantee your stored `AnnData` object will be fully readable.
And even if it is fully readable, it will almost certainly be much slower to read.

There are two ways of opening remote `zarr` stores from the `zarr-python` package, {class}`zarr.storage.FsspecStore` and {class}`zarr.storage.ObjectStore`, and both can be used with `read_lazy`.
[`obstore` claims] to be more performant out-of-the-box, but notes that this claim has not been benchmarked with the `uvloop` event loop, which itself claims to be 2× more performant than the default event loop for `python`.

## Local data

Local data generally poses a different set of challenges.
First, write speeds can be somewhat slow and second, the creation of many small files on a file system can slow down a filesystem.
For the "many small files" problem, `zarr` has introduced {ref}`sharding <zarr:user-guide-sharding>` in the v3 file format.
Sharding requires knowledge of the array element you are writing (such as shape or data type), though, and therefore you will need to use {func}`anndata.experimental.write_dispatched` to use sharding.
For example, you cannot shard a 1D array with `shard` sizes `(256, 256)`.
Here is a short example, although you should tune the sizes to your own use-case and also use the compression that makes the most sense for you:

```python
import zarr
import anndata as ad
from collections.abc import Mapping
from typing import Any

ad.settings.zarr_write_format = 3 # Absolutely crucial! Sharding is only for the v3 file format!

def write_sharded(group: zarr.Group, adata: ad.AnnData):
    def callback(
        func: ad.experimental.Write,
        g: zarr.Group,
        k: str,
        elem: ad.typing.RWAble,
        dataset_kwargs: Mapping[str, Any],
        iospec: ad.experimental.IOSpec,
    ):
        if iospec.encoding_type in {"array"}:
            dataset_kwargs = {
                "shards": tuple(int(2 ** (16 / len(elem.shape))) for _ in elem.shape),
                **dataset_kwargs,
            }
            dataset_kwargs["chunks"] = tuple(i // 2 for i in dataset_kwargs["shards"])
        elif iospec.encoding_type in {"csr_matrix", "csc_matrix"}:
            dataset_kwargs = {"shards": (2**16,), "chunks": (2**8,), **dataset_kwargs}
        func(g, k, elem, dataset_kwargs=dataset_kwargs)

    return ad.experimental.write_dispatched(group, "/", adata, callback=callback)
```

However, `zarr-python` can be slow with sharding throughput as well as writing throughput.
Thus if you wish to speed up either writing, sharding, or both (or receive a modest speed-boost for reading), a bridge to the `zarr` implementation in Rust {doc}`zarrs-python <zarrs:index>` can help with that (see the [zarr-benchmarks]):

```
uv pip install zarrs
```

```python
import zarr
import zarrs
zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})
```

However, this pipeline is not compatible with all types of zarr store, especially remote stores and there are limitations on where rust can give a performance boost for indexing.
We therefore recommend this pipeline for writing full datasets and reading contiguous regions of said written data.

## Codecs

The default `zarr-python` v3 codec for the v3 format is no longer `blosc` but `zstd`.
While `zstd` is more widespread, you may find its performance to not meet your old expectations.
Therefore, we recommend passing in the {class}`zarr.codecs.BloscCodec` to `compressor` on {func}`~anndata.AnnData.write_zarr` if you wish to return to the old behavior.

There is currently a bug with `numcodecs` that prevents data written from other non-numcodecs `zstd` implementations from being read in by the default zarr pipeline (to which the above rust pipeline falls back if it cannot handle a datatype or indexing scheme, like `vlen-string`): {issue}`zarr-developers/numcodecs#424`.
Thus is may be advisable to use `BloscCodec` with `zarr` v3 file format data if you wish to use the rust-accelerated pipeline until this issue is resolved.

The same issue with `zstd` applies to data that may eventually be written by the GPU `zstd` implementation (see below).

## Dask

Zarr v3 should be compatible with dask, although the default behavior is to use zarr's chunking for dask's own. With sharding, this behavior may be undesirable as shards can often contain many small chunks, thereby slowing down i/o as dask will need to index into the zarr store for every chunk.  Therefore it may be better to customize this behavior by passing `chunks=my_zarr_array.shards` as an argument to {func}`dask.array.from_zarr` or similar.

## GPU i/o

At the moment, it is unlikely your `anndata` i/o will work if you use {ref}`zarr.config.enable_gpu <zarr:user-guide-gpu>`.
It's *possible* dense data i/o i.e., using {func}`anndata.io.read_elem` will work as expected, but this functionality is untested – sparse data, awkward arrays, and dataframes will not.
`kvikio` currently provides a {class}`kvikio.zarr.GDSStore` although there are no working compressors at the moment exported from the `zarr-python` package (work is underway for `Zstd`: {pr}`zarr-developers/zarr-python#2863`.

We anticipate enabling officially supporting this functionality officially for dense data, sparse data, and possibly awkward arrays in the next minor release, 0.13.

## Asynchronous i/o

At the moment, `anndata` exports no `async` functions.
However, `zarr-python` has a fully `async` API and provides its own event-loop so that users like `anndata` can interact with a synchronous API while still beenfitting from `zarr-python`'s asynchronous functionality under that API.
We anticipate providing `async` versions of {func}`anndata.io.read_elem` and {func}`anndata.experimental.read_dispatched` so that users can download data asynchronously without using the `zarr-python` event loop.
We also would like to create an asynchronous partial reader to enable iterative streaming of a dataset.

[`obstore` claims]: https://developmentseed.org/obstore/latest/performance
[zarr-benchmarks]: https://github.com/LDeakin/zarr_benchmarks
