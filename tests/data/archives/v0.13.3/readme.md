These files were written with

```bash
uvx '--with=anndata[dask]==0.13.3.post0' '--with=zarr==3.3' '--with=pytest==9.1.1' '--with=awkward==2.13.0' python -c '
import zarr
import h5py
from anndata.tests.helpers import gen_adata
import anndata as ad
ad.settings.zarr_write_format = 2
adata = gen_adata((10, 20))
adata.write_zarr(zarr.storage.ZipStore("tests/data/archives/v0.13.3/adata.zarr.zip", mode="w"))
adata.write_h5ad("tests/data/archives/v0.13.3/adata.h5ad")'
```a
