"""\
Utility functions for AnnData._repr_html_()
"""
import uuid
from functools import lru_cache, singledispatch
from html import escape
import contextlib
from typing import Union, Literal, TypedDict, Mapping  # Generic ABCs
from importlib.resources import read_binary
from scipy import sparse
import pandas as pd
import numpy as np


STATIC_FILES = (
    ("anndata.static.html", "icons-svg-inline.html"),
    ("anndata.static.css", "style.css"),
)

Options = Literal[
    "display_max_rows",
    "display_values_threshold",
    "display_style",
    "display_width",
    "display_expand_attrs",
    "display_expand_coords",
    "display_expand_data_vars",
    "display_expand_data",
]


class T_Options(TypedDict):
    display_max_rows: int
    display_values_threshold: int
    display_style: Literal["text", "html"]
    display_width: int
    display_expand_attrs: Literal["default", True, False]
    display_expand_coords: Literal["default", True, False]
    display_expand_data_vars: Literal["default", True, False]
    display_expand_data: Literal["default", True, False]


OPTIONS: T_Options = {
    "display_max_rows": 12,
    "display_values_threshold": 200,
    "display_style": "html",
    "display_width": 100,
    "display_expand_attrs": "default",
    "display_expand_coords": "default",
    "display_expand_data_vars": "default",
    "display_expand_data": "default",
}


@lru_cache(None)
def _load_static_files():
    """Lazily load the resource files into memory the first time they are needed"""
    return [
        read_binary(package, resource).decode("utf-8")
        for package, resource in STATIC_FILES
    ]


@singledispatch
def _html_format(x):
    return escape(str(x))


@_html_format.register(pd.DataFrame)
def _html_format_df(x: pd.DataFrame):
    return x._repr_html_()


@_html_format.register(np.ndarray)
def _html_format_np(x: np.ndarray):
    # default to lower precision so a full (abbreviated) line can fit on
    # one line with the default display_width
    options = {
        "precision": 6,
        "linewidth": OPTIONS["display_width"],
        "threshold": OPTIONS["display_values_threshold"],
    }
    if x.ndim < 3:
        edgeitems = 3
    elif x.ndim == 3:
        edgeitems = 2
    else:
        edgeitems = 1
    options["edgeitems"] = edgeitems
    with set_numpy_options(**options):
        return repr(x)


@_html_format.register(sparse.spmatrix)
def _html_format_sp(x: sparse.spmatrix):
    return f"<pre>{escape(x.__repr__())}</pre>"


@singledispatch
def _dtype_repr(x):
    return escape(str(x.dtype))


@_dtype_repr.register(pd.DataFrame)
def _dtype_repr_pd(x: pd.DataFrame):
    return f"DataFrame {escape(str(x.shape))}"


# TODO: document and mention
# https://github.com/pydata/xarray/blob/8e9a9fb390f8f0d27a017a7affd8d308d2317959/xarray/core/formatting.py#L32
def maybe_truncate(obj, max_width=500):
    s = str(obj)
    if len(s) > max_width:
        s = s[: (max_width - 3)] + "..."
    return s


@singledispatch
def _inline_format(x, max_width):
    return maybe_truncate(escape(str(x)), max_width=max_width)


@contextlib.contextmanager
def set_numpy_options(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally:
        np.set_printoptions(**original)


def _obj_repr(obj, header_components, sections):
    """Return HTML repr of an anndata object.

    If CSS is not injected (untrusted notebook), fallback to the plain text repr.

    """
    header = f"<div class='ad-header'>{''.join(h for h in header_components)}</div>"
    sections = "".join(f"<li class='ad-section-item'>{s}</li>" for s in sections)

    icons_svg, css_style = _load_static_files()
    return (
        "<div>"
        f"{icons_svg}<style>{css_style}</style>"
        f"<pre class='ad-text-repr-fallback'>{escape(repr(obj))}</pre>"
        "<div class='ad-wrap' style='display:none'>"
        f"{header}"
        f"<ul class='ad-sections'>{sections}</ul>"
        "</div>"
        "</div>"
    )


def _summarize_attrs(attrs):
    attrs_dl = "".join(
        f"<dt><span>{escape(str(k))} :</span></dt>" f"<dd>{escape(str(v))}</dd>"
        for k, v in attrs.items()
    )

    return f"<dl class='ad-attrs'>{attrs_dl}</dl>"


def _icon(icon_name):
    # TODO: Change names etc.
    # icon_name should be defined in anndata/static/html/icon-svg-inline.html
    return (
        "<svg class='icon ad-{0}'>"
        "<use xlink:href='#{0}'>"
        "</use>"
        "</svg>".format(icon_name)
    )


def _summarize_item_html(
    name: str, x: Union[pd.DataFrame, np.ndarray, sparse.spmatrix], attrs=None
):
    """
    Summarizes x
    """
    # TODO: more detail

    # TODO: learn what this is
    cssclass_idx = ""

    dims_str = f"({', '.join(escape(str(dim)) for dim in x.shape)})"
    name = escape(str(name))
    dtype = _dtype_repr(x)

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())
    disabled = "" if attrs else "disabled"

    preview = _inline_format(x, 35)
    data_repr = _html_format(x)
    attrs_ul = _summarize_attrs(attrs if attrs else {})

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    return (
        f"<div class='ad-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<div class='ad-var-dims'>{dims_str}</div>"
        f"<div class='ad-var-dtype'>{dtype}</div>"
        f"<div class='ad-var-preview ad-preview'>{preview}</div>"
        f"<input id='{attrs_id}' class='ad-var-attrs-in' "
        f"type='checkbox' {disabled}>"
        f"<label for='{attrs_id}' title='Show/Hide attributes'>"
        f"{attrs_icon}</label>"
        f"<input id='{data_id}' class='ad-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>"
        f"{data_icon}</label>"
        f"<div class='ad-var-attrs'>{attrs_ul}</div>"
        f"<div class='ad-var-data'>{data_repr}</div>"
    )


def _collapsible_section(
    name,
    inline_details="",
    details="",
    n_items=None,
    enabled=True,
    collapsed=False,
    **kwargs,
):
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())

    has_items = n_items is not None and n_items
    n_items_span = "" if n_items is None else f" <span>({n_items})</span>"
    enabled = "" if enabled and has_items else "disabled"
    collapsed = "" if collapsed or not has_items else "checked"
    tip = " title='Expand/collapse section'" if enabled else ""

    return (
        f"<input id='{data_id}' class='ad-section-summary-in' "
        f"type='checkbox' {enabled} {collapsed}>"
        f"<label for='{data_id}' class='ad-section-summary' {tip}>"
        f"{name}:{n_items_span}</label>"
        f"<div class='ad-section-inline-details'>{inline_details}</div>"
        f"<div class='ad-section-details'>{details}</div>"
    )


def _summarize_mapping_html(x: Mapping):
    vars_li = "".join(
        f"<li class='ad-var-item'>{_summarize_item_html(k, v)}</li>"
        for k, v in x.items()
    )
    return f"<ul class='ad-var-list'>{vars_li}</ul>"


def _create_sections_from_conf(ad_obj, sections_conf):
    sections = []
    for k, v in sections_conf.items():
        # TODO: maybe use enums?
        if v["section_type"] == "single":
            sections.append(
                _collapsible_section(
                    name=k,
                    details=_summarize_item_html(
                        name=k, x=getattr(ad_obj, k), attrs=v.get("attrs", {})
                    ),
                    **v["args"],
                )
            )
        elif v["section_type"] == "mapping":
            sections.append(
                _collapsible_section(
                    name=k,
                    details=_summarize_mapping_html(x=getattr(ad_obj, k)),
                    **v["args"],
                )
            )
        else:
            raise NotImplementedError()

    return sections


def _create_anndata_repr(ad_obj: "AnnData"):
    obj_type = f"anndata.{type(ad_obj).__name__}"

    header_components = [f"<div class='ad-obj-type'>{escape(obj_type)}</div>"]

    sections_conf = {
        # sections consisting of single items like matrices
        "X": {
            "section_type": "single",
            "attrs": {
                "Is backed": "Nowhere"
                if not ad_obj.isbacked
                else escape(ad_obj.file.filename)
            },
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": ad_obj.n_obs,
            },
        },
        "obs": {
            "section_type": "single",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": ad_obj.n_obs,
            },
        },
        "var": {
            "section_type": "single",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": ad_obj.n_vars,
            },
        },
        # dict like sections
        "obsm": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": len(ad_obj.obsm.keys()),  # TODO: Handle when obsm is ndarray
            },
        },
        "varm": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": len(ad_obj.varm.keys()),
            },
        },
        "obsp": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": len(ad_obj.obsp.keys()),
            },
        },
        "varp": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": len(ad_obj.varp.keys()),
            },
        },
    }
    if ad_obj.layers:
        sections_conf["layers"] = {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": False,
                "n_items": len(ad_obj.layers.keys()),
            },
        }

    sections = _create_sections_from_conf(ad_obj, sections_conf)
    return _obj_repr(ad_obj, header_components, sections)
