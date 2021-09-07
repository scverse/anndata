import pandas as pd
from typing import Optional, Union, List, Dict

from .anndata import AnnData

GroupsList = List[str]
GroupsLists = List[GroupsList]
GroupsDict = Dict[str, GroupsList]


def split(
    adata: AnnData,
    by: Union[str, pd.Series, list],
    groups: Optional[Union[GroupsList, GroupsLists, GroupsDict]] = None,
    others_key: Optional[str] = None,
    axis: int = 0,
    copy: bool = False,
) -> Dict[str, AnnData]:
    """\
    Split adata by obs key.
    Params
    ------
    adata
        AnnData object to split.
    by
        A key from `.obs` or `.var` depending on `axis`,
        pandas Series or a list with groups to use for splitting `adata`.
    groups
        Specifies which groups to select from adata `by`.
        If `None`, `adata` will be split with all values from `by`.
        It can be a list of groups' names, a list of lists of groups to aggregate,
        a dict of lists of groups to aggregate.
    others_key
        If not `None`, the returned dict will have an additional key with
        this name and the adata subset object with all groups not specified
        in `groups` as a value.
    axis
        Axis of `adata` to split by.
    copy
        If `True`, all split AnnData objects are copied; otherwise,
        the returned dict will have views.
    Returns
    -------
    A dictionay with AnnData objects.
    If names of the split AnnData objects are not specified
    (a list of lists in `groups`), they will be created by joining the groups' names.
    Examples
    --------
    >>> import scanpy as sc
    >>> import pandas as pd
    >>> from anndata import split
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> # Split by all values in an `.obs` key:
    >>> adatas = split(adata, 'bulk_labels')
    >>> adatas #doctest: +ELLIPSIS
    {'CD14+ Monocyte': View of AnnData object with n_obs × n_vars = 129 × 765...}
    >>> # Select only specific groups from `.obs` key:
    >>> adatas = split(adata, 'bulk_labels', ['CD14+ Monocyte', 'CD34+'])
    >>> adatas #doctest: +ELLIPSIS
    {'CD14+ Monocyte': View of AnnData object with n_obs × n_vars = 129 × 765...}
    >>> # Aggreagte some groups from `.obs` key, put all others to `others_key`:
    >>> adatas = split(adata, 'bulk_labels', dict(some=['CD14+ Monocyte', 'CD34+']),
    ...                others_key='others')
    >>> adatas #doctest: +ELLIPSIS
    {'some': View of AnnData object with n_obs × n_vars = 142 × 765...}
    >>> # Split by axis 1 (var axis) passing the Series from the binned continuous variable
    >>> var_cats = pd.cut(adata.var.means, 4).cat.rename_categories(str)
    >>> adatas = split(adata, var_cats, axis=1)
    >>> adatas #doctest: +ELLIPSIS
    {'(0.0233, 1.135]': View of AnnData object with n_obs × n_vars = 700 × 535...}
    """
    if axis not in (0, 1):
        raise ValueError("axis should be 0 or 1 only.")

    attr_axis, attr_str = (adata.obs, ".obs") if axis == 0 else (adata.var, ".var")

    if isinstance(by, str):
        if by not in attr_axis:
            raise ValueError(f"No {by} in {attr_str}.")
        attr_by = attr_axis[by]
    elif isinstance(by, list):
        attr_by = pd.Series(by)
    elif isinstance(by, pd.Series):
        attr_by = by
    else:
        raise ValueError("by should be a string, pandas Series or a list")

    if len(attr_by) != len(attr_axis.index):
        raise ValueError(f"by should have the same length as {attr_str}.")

    select = [slice(None), slice(None)]

    adatas = {}

    if groups is None:
        groups = attr_by.unique()

    groups_dict = {}
    all_values = []
    for group in groups:
        if isinstance(groups, dict):
            values = groups[group]
        else:
            values = group

        if isinstance(group, list):
            name = "-".join(str(e) for e in group)
        else:
            name = group

        values = values if isinstance(values, list) else [values]

        groups_dict[name] = values
        all_values += values

    use_others_key = others_key is not None
    # need to create dict before checking that others_key
    # is not among the passed groups
    if use_others_key and others_key in groups_dict:
        raise ValueError(
            f"others_key={others_key} coincides with a key in the passed groups."
        )

    for group, values in groups_dict.items():
        select[axis] = attr_by.isin(values)
        idx = tuple(select)
        adatas[group] = adata[idx].copy() if copy else adata[idx]

    # should be last
    if use_others_key:
        mask = ~attr_by.isin(all_values)
        if sum(mask) > 0:
            select[axis] = mask
            idx = tuple(select)
            adatas[others_key] = adata[idx].copy() if copy else adata[idx]

    return adatas
