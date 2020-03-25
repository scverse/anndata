from functools import singledispatch
from collections.abc import Mapping

from ipytree import Tree, Node

import numpy as np
from scipy import sparse
import pandas as pd

import anndata as ad


def anndata_descr(x: ad.AnnData) -> str:
    descr = f"AnnData object with {x.n_obs} obs × {x.n_vars} vars"
    if x.is_view:
        descr = "View of " + descr
    if x.isbacked:
        descr += f" backed at '{x.filename}'"
    return descr


def format_name(name):
    return f"<code><b>{name}</code></b>"


def as_tree(x: ad.AnnData, **kwargs) -> Tree:
    tree = Tree(stripes=True)
    root = Node(anndata_descr(x), icon="sitemap")
    tree.add_node(root)
    elnames = ["X", "obs", "var", "layers", "obsm", "varm", "obsp", "varp", "uns"]
    for elname in elnames:
        el = getattr(x, elname)
        if isinstance(el, Mapping) and len(el) == 0:
            continue
        root.add_node(as_node(el, elname, **kwargs))
    return tree


@singledispatch
def as_node(x, name: str, **kwargs):
    raise NotImplementedError(f"Not implemented for type {type(x)}")


@as_node.register
def as_node_mapping(x: Mapping, name: str, **kwargs):
    node = Node(format_name(name), icon="sitemap", **kwargs)
    for k, v in x.items():
        node.add_node(as_node(v, name=k, **kwargs))
    return node


@as_node.register
def as_node_alignedmapping(
    x: ad._core.aligned_mapping.AlignedMapping, name: str, **kwargs
):
    node = Node(
        f"{format_name(name)}: <code>Mapping</code> aligned to <code>{x.dims}</code> with <code>{len(x)}</code> entries.",
        icon="clone",
        **kwargs,
    )
    for k, v in x.items():
        node.add_node(as_node(v, name=k, **kwargs))
    return node


# @as_node.registed(ad.SparseDataset)
# @as_node.register(h5py.DataSet)
@as_node.register(np.ndarray)
def as_node_arr(x: np.ndarray, name: str, icon="square", **kwargs):
    return Node(
        f"{format_name(name)}: <code>{type(x).__name__}</code> of shape <code>{'×'.join(map(str, x.shape))}</code> and dtype <code>{x.dtype}</code>",
        icon=icon,
        **kwargs,
    )


@as_node.register
def as_node_df(x: pd.DataFrame, name: str, **kwargs):
    node = Node(
        f"{format_name(name)}: <code>{type(x).__name__}</code> of shape <code>{'×'.join(map(str, x.shape))}</code>",
        icon="table",
        **kwargs,
    )
    for k, v in x.items():
        node.add_node(as_node(v, name=k, **kwargs))
    return node


@as_node.register
def as_node_series(x: pd.Series, name: str, **kwargs):
    return Node(
        f"{format_name(name)}: <code>{type(x).__name__}</code> with dtype <code>{x.dtype}</code>",
        icon="columns",
        **kwargs,
    )


@as_node.register
def as_node_sparse(x: sparse.spmatrix, name: str, **kwargs):
    return as_node_arr(x, name, **kwargs)


@as_node.register(type(None))
@as_node.register(int)
@as_node.register(float)
@as_node.register(bool)
@as_node.register(np.bool_)
@as_node.register(np.number)
@as_node.register(np.string_)
@as_node.register(np.str_)
@as_node.register(np.str)
def as_node_scalar(x, name: str, **kwargs):
    return Node(
        f"{format_name(name)}: <code>{repr(x)}</code> (<code>{type(x).__name__}</code>)",
        icon="minus",
        **kwargs,
    )
