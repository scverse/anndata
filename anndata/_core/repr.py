#
# Utility functions for AnnData._repr_html_()
#

from typing import Tuple, Iterable
from numbers import Integral, Real, Complex
from warnings import warn
import numpy as np
import pandas as pd

def repr_html(adata: "AnnData", expand=None) -> str:
    """
    HTML formatter for the AnnData object.

    This formatter has an optional argument `expand`
    which is a 2-bit flag:
    10 - expand AnnData slots
    01 - expand individul slots
    """

    if expand is None:
        expand = 0b10

    data = ""

    data += "<div class='block-anndata'><div>"
    data += "<details{}>".format(" open" if (expand & 0b10) >> 1 else "")
    data += "<summary class='summary-anndata'><div class='title title-anndata'>{}</div><span class='hl-dim'>{} &times {}</span></summary>".format(
        "AnnData", *adata.shape
    )

    # General object properties
    data += (
        "<span>{} object <span class='hl-dim'>{} obs &times; {} var</span></span>".format(
            type(adata).__name__, *(adata.shape)
        )
    )
    if adata.isbacked:
        data += "<br>&#8627; <span>backed at <span class='hl-file'>{}</span></span>".format(
            adata.file.filename
        )

    data += "<br>"

    # X
    data += block_matrix(adata, "X", "Matrix")
    # Layers
    data += details_block_table(adata, "layers", "Layers", expand & 0b01, dims=False)
    # Metadata
    data += details_block_table(adata, "obs", "Metadata", expand & 0b01)
    # Embeddings
    data += details_block_table(adata, "obsm", "Embeddings", expand & 0b01)
    # Distances
    data += details_block_table(adata, "obsp", "Distances", expand & 0b01, square=True)
    # Miscellaneous (unstructured)
    data += details_block_table(adata, "uns", "Miscellaneous", expand & 0b01)

    data += "</details>"
    data += "</div></div>"

    data += "<br/>"
    full = "".join((ANNDATA_CSS, data))
    return full


def maybe_module_class(obj, sep=".", builtins=False) -> Tuple[str, str]:
    m, cl = "", obj.__class__.__name__
    try:
        m += obj.__class__.__module__
        if m == "builtins" and not builtins:
            m = ""
        else:
            m += sep
    except:
        m += ""
    return (m, cl)


def format_values(x):
    s = ""
    if not isinstance(x, Iterable):
        s += f"{x}"
    elif isinstance(x, pd.DataFrame):
        s += "DataFrame ({} x {})".format(*x.shape)
    elif hasattr(x, "keys") and hasattr(x, "values") and not hasattr(x, "shape"):
        ks = ",".join(x.keys())
        if "," not in ks:  # only 1 element
            vs = ",".join(map(format_values, x.values()))
            s += ks + ": " + vs
        else:
            s += ks
    elif isinstance(x, str):
        s += x
    else:
        x = x[: min(100, len(x))]
        if hasattr(x, "shape"):
            if isinstance(x, np.ndarray):
                x = x.flat
            elif isinstance(x, pd.Series):
                x = x.to_numpy()
            else:
                warn(f"got unknown array type {type(x)}, don't know how handle it.")
                return type(x)
        if isinstance(x[0], Integral):
            s += ",".join([f"{i}" for i in x])
        elif isinstance(x[0], Real):
            s += ",".join([f"{i:.2f}" for i in x])
        elif isinstance(x[0], Complex):
            warn("got complex number, don't know how to handle it")
        elif isinstance(x[0], Iterable):
            s += ",".join(map(format_values, x))
        s = s[:50]
        while s[-1] != ",":
            s = s[:-1]
        s += "..."
    return s


def block_matrix(data, attr, name):
    obj = getattr(data, attr)
    s = ""
    s += "<div class='title title-attr'>{}</div><span class='hl-dim'>.{}</span>".format(name, attr)
    s += "<div>"
    s += """
            <span class='hl-types'>{}</span> <span>&nbsp;&nbsp;&nbsp;<span class='hl-import'>{}</span>{}</span>
         """.format(
        obj.dtype, *maybe_module_class(obj)
    )
    s += "</table></div>"
    return s


def details_block_table(data, attr, name, expand=0, dims=True, square=False):
    obj = getattr(data, attr)
    s = ""
    # DataFrame
    if isinstance(obj, pd.DataFrame):
        s += "<details{}>".format(" open" if expand else "")
        s += "<summary><div class='title title-attr'>{}</div><span class='hl-dim'>.{}</span><span class='hl-size'>{} element{}</span></summary>".format(
            name, attr, obj.shape[1], "s" if obj.shape[1] != 1 else ""
        )
        if obj.shape[1] > 0:
            s += "<div><table>"
            s += "\n".join(
                [
                    """<tr>
                                <td class='col-index'>{}</td>  <td class='hl-types'>{}</td>  <td class='hl-values'>{}</td>
                            </tr>""".format(
                        attr_key, obj[attr_key].dtype, format_values(obj[attr_key])
                    )
                    for attr_key in obj.columns
                ]
            )
            s += "</table></div>"
        else:
            s += f"<span class='hl-empty'>No {name.lower()}</span>"
        s += "</details>"
    # Dict-like object
    elif hasattr(obj, "keys") and hasattr(obj, "values") and name != "Miscellaneous":
        s += "<details{}>".format(" open" if expand else "")
        s += "<summary><div class='title title-attr'>{}</div><span class='hl-dim'>.{}</span><span class='hl-size'>{} element{}</span></summary>".format(
            name, attr, len(obj), "s" if len(obj) != 1 else ""
        )
        if len(obj) > 0:
            s += "<div><table>"
            if square:  # e.g. distance matrices in .obsp
                s += "\n".join(
                    [
                        """<tr>
                                       <td class='col-index'>{}</td>  <td class='hl-types'>{}</td>  <td><span class='hl-import'>{}</span>{}</td>
                                   </tr>""".format(
                            attr_key, obj[attr_key].dtype, *maybe_module_class(obj[attr_key])
                        )
                        for attr_key in obj.keys()
                    ]
                )
            else:  # e.g. embeddings in .obsm
                s += "\n".join(
                    [
                        """<tr>
                                       <td class='col-index'>{}</td>  <td class='hl-types'>{}</td>  <td><span class='hl-import'>{}</span>{}</td>  <td class='hl-dims'>{}</td>
                                   </tr>""".format(
                            attr_key,
                            obj[attr_key].dtype,
                            *maybe_module_class(obj[attr_key]),
                            f"{obj[attr_key].shape[1]} dims"
                            if len(obj[attr_key].shape) > 1 and dims
                            else "",
                        )
                        for attr_key in obj.keys()
                    ]
                )

            s += "</table></div>"
        else:
            s += f"<span class='hl-empty'>No {name.lower()}</span>"
        s += "</details>"
    elif hasattr(obj, "file"):  # HDF5 dataset
        s += "<details{}>".format(" open" if expand else "")
        s += "<summary><div class='title title-attr'>{}</div><span class='hl-dim'>.{}</span><span class='hl-size'>{} elements</span></summary>".format(
            name, attr, len(obj)
        )
        s += "<div><table>"
        s += """<tr>
                <td class='hl-types'>{}</td>  <td><span class='hl-import'>{}</span>{}</td>
                </tr>""".format(
            obj.dtype, *maybe_module_class(obj)
        )
        s += "</table></div>"
        s += "</details>"
    else:  # Unstructured
        s += "<details{}>".format(" open" if expand else "")
        s += "<summary><div class='title title-attr'>{}</div><span class='hl-dim'>.{}</span><span class='hl-size'>{} elements</span></summary>".format(
            name, attr, len(obj)
        )
        if len(obj) > 0:
            s += "<div><table>"
            s += "\n".join(
                [
                    """<tr>
                                <td class='col-index'>{}</td>  <td><span class='hl-import'>{}</span>{}</td>  <td class='hl-dims'>{} element{}</td>  <td class='hl-values'>{}</td>
                            </tr>""".format(
                        attr_key,
                        *maybe_module_class(obj[attr_key]),
                        len(obj[attr_key]),
                        "s" if len(obj[attr_key]) != 1 else "",
                        format_values(obj[attr_key]),
                    )
                    for attr_key in obj.keys()
                ]
            )
            s += "</table></div>"
        else:
            s += f"<span class='hl-empty'>No {name.lower()}</span>"
        s += "</details>"
    return s


ANNDATA_CSS = """<style>
.hl-dim, .hl-size, .hl-values, .hl-types, .hl-dims {
  color: #777777;
}
.hl-dim::before, .hl-size::before {
  content: "\\00a0\\00a0\\00a0";
}
.hl-values {
  font-family: monospace;
}
.hl-file {
  background-color: #EEEEEE;
  border-radius: .5rem;
  padding: .2rem .4rem;
  color: #555555;
}
.hl-empty {
  color: #999999;
}
.hl-import {
  color: #777777;
}
.block-anndata {
  display: block;
  margin: 0 2rem;
}
.block-anndata .title {
  display: inline-block;
  font-weight: 600;
  color: #555555;
}
.block-anndata .title-anndata {
  font-size: 1.2rem;
  color: #ee733c;
  padding: 0 .5rem;
}
.block-anndata .title-attr {
  font-size: 1.0rem;
  padding-top: .2rem;
}
.block-anndata summary {
  cursor: pointer;
  list-style: none;
}
.block-anndata summary::-webkit-details-marker {
  display: none;
}
.block-anndata details > summary::before {
  content: '\u2295';
}
.block-anndata details[open] > summary::before {
  content: '\u2296';
}
.block-anndata table tr {
  background-color: transparent !important;
}
.block-anndata table tr:hover {
  background-color: #ee733c55 !important;
}
.block-anndata .col-index {
  text-align: left !important;
}
.block-anndata .summary-anndata {
  margin-left: -2rem;
}
.block-anndata > .summary-anndata:hover {
  background-color: #ee733c55;
}
.block-anndata .summary-anndata::before {
  color: #ee733c;
  content: '\u25cf';
}
</style>"""
