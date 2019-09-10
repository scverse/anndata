import warnings

import numpy as np
import pandas as pd


def version(package):
    try:
        from importlib.metadata import version
    except ImportError:
        from importlib_metadata import version
    return version(package)


def warn_flatten():
    warnings.warn(
        "In anndata v0.7+, arrays contained within an AnnData object will "
        "maintain their dimensionality. For example, prior to v0.7 `adata[0, 0].X`"
        " returned a scalar and `adata[0, :]` returned a 1d array, post v0.7 they"
        " will return two dimensional arrays. If you would like to get a one "
        "dimensional array from your AnnData object, consider using the "
        "`adata.obs_vector`, `adata.var_vector` methods or accessing the array"
        " directly.",
        FutureWarning,
        stacklevel=2,
    )


def _from_fixed_length_strings(value):
    """Convert from fixed length strings to unicode.

    For backwards compatability with older h5ad and zarr files.
    """
    new_dtype = []
    for dt in value.dtype.descr:
        dt_list = list(dt)
        dt_type = dt[1]
        if isinstance(dt_type, tuple):  # vlen strings, could probably match better
            dt_list[1] = "O"
            new_dtype.append(tuple(dt_list))
        elif issubclass(np.dtype(dt_type).type, np.string_):
            dt_list[1] = 'U{}'.format(int(dt_type[2:]))
            new_dtype.append(tuple(dt_list))
        else:
            new_dtype.append(dt)
    return value.astype(new_dtype)


def _clean_uns(d: dict):
    """Compat function for when categorical keys were stored in uns."""
    k_to_delete = []
    for k, v in d.get("uns", {}).items():
        if k.endswith('_categories'):
            k_stripped = k.replace('_categories', '')
            if isinstance(v, (str, int)):  # fix categories with a single category
                v = [v]
            for ann in ['obs', 'var']:
                if k_stripped in d[ann]:
                    d[ann][k_stripped] = pd.Categorical.from_codes(
                        codes=d[ann][k_stripped].values, categories=v
                    )
                    k_to_delete.append(k)
    for k in k_to_delete:
        del d["uns"][k]
