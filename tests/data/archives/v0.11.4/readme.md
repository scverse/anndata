These files were written with

```bash
uvx '--with=anndata==0.11.4' '--with=zarr<3' python -c '
import zarr
from anndata import AnnData
adata = AnnData(shape=(10, 20))
adata.write_zarr(zarr.ZipStore("tests/data/archives/v0.11.4/adata.zarr.zip"))
adata.write_h5ad("tests/data/archives/v0.11.4/adata.h5ad")'
```
