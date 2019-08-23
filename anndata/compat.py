import warnings


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
        stacklevel=2
    )
