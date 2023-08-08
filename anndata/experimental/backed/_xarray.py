import xarray as xr
from anndata._core.index import Index, _subset
from anndata._core.views import as_view


def get_index_dim(ds):
    assert (
        len(ds.dims) == 1
    ), f"xarray Dataset should not have more than 1 dims, found {len(ds)}"
    return list(ds.dims.keys())[0]


class Dataset2D(xr.Dataset):
    @property
    def shape(
        self,
    ):  # aligned mapping classes look for this for DataFrames so this ensures usability with e.g., obsm
        return [self.dims[get_index_dim(self)], len(self)]


@_subset.register(Dataset2D)
def _(a: xr.DataArray, subset_idx: Index):
    key = get_index_dim(a)
    if (
        isinstance(subset_idx, tuple) and len(subset_idx) == 1
    ):  # xarray seems to have some code looking for a second entry in tuples
        return a.isel(**{key: subset_idx[0]})
    return a.isel(**{key: subset_idx})


@as_view.register(Dataset2D)
def _(a: Dataset2D, view_args):
    return a
