"""\
Utility functions for AnnData._repr_html_()
TODO: License stuff for Xarray codes here.
"""
import sys
import uuid
from functools import lru_cache, singledispatch
from html import escape
from typing import Union, Mapping
from importlib.resources import read_binary
from scipy import sparse
import pandas as pd
import numpy as np
from . import anndata

if sys.version_info >= (3, 8):
    from typing import TypedDict, Literal
else:
    from typing_extensions import TypedDict, Literal


STATIC_FILES = (
    ("anndata.static.html", "icons-svg-inline.html"),
    ("anndata.static.css", "style.css"),
)

Options = Literal[
    "display_max_rows",
    "display_values_threshold",
    "display_style",
    "display_width",
    "display_expand_attrs",  # the option with the page icon
    "display_expand_data",  # the option with the database icon
    "display_expand_mapping_section",  # the sections where data is mapping type e.g. obsm, varm, obsp...
    "display_expand_single_item_section",  # the sections where data is not mapping type e.g. X, var, obs...
]


class T_Options(TypedDict):
    display_max_rows: int
    display_values_threshold: int
    display_style: Literal["text", "html"]
    display_width: int
    display_expand_attrs: Literal["default", True, False]
    display_expand_data: Literal["default", True, False]
    display_expand_mapping_section: Literal["default", True, False]
    display_expand_single_item_section: Literal["default", True, False]


OPTIONS: T_Options = {
    "display_max_rows": 12,
    "display_values_threshold": 50,
    "display_style": "html",
    "display_width": 75,
    "display_expand_attrs": "default",
    "display_expand_data": "default",
    "display_expand_mapping_section": "default",
    "display_expand_single_item_section": "default",
}


def _get_boolean_with_default(option: Options, default: bool) -> bool:
    global_choice = OPTIONS[option]

    if global_choice == "default":
        return default
    elif isinstance(global_choice, bool):
        return global_choice
    else:
        raise ValueError(
            f"The global option {option} must be one of True, False or 'default'."
        )


@lru_cache(None)
def _load_static_files():
    """Lazily load the resource files into memory the first time they are needed"""
    return [
        read_binary(package, resource).decode("utf-8")
        for package, resource in STATIC_FILES
    ]


# TODO: document and mention
# https://github.com/pydata/xarray/blob/8e9a9fb390f8f0d27a017a7affd8d308d2317959/xarray/core/formatting.py#L32
def _maybe_truncate(obj, max_width=500):
    s = str(obj)
    if len(s) > max_width:
        s = s[: (max_width - 3)] + "..."
    return s


"""
HTML FORMAT FUNCTIONS
These functions give an representation to be shown in the html summary.
"""


@singledispatch
def _html_format(x):
    return f"<pre>{escape(str(x))}</pre>"


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
    with np.printoptions(**options):
        return f"<pre>{escape(x.__repr__())}</pre>"


@_html_format.register(sparse.spmatrix)
def _html_format_sp(x: sparse.spmatrix):
    return f"<pre>{escape(x.__repr__())}</pre>"


"""
Dimension repr FUNCTIONS
Gives the shape of different data-structures based on their type.
"""



@singledispatch
def _dim_repr(x):
    if hasattr(x, "shape"):
        return _dim_repr_shape(x)
    return "Unstructured"


@_dim_repr.register(pd.DataFrame)
@_dim_repr.register(sparse.spmatrix)
@_dim_repr.register(np.ndarray)
def _dim_repr_shape(x):
    return f"({', '.join(escape(str(dim)) for dim in x.shape)})"


"""
Type repr FUNCTIONS
Gives a short representation of the type of the data structure.
"""


@singledispatch
def _type_repr(x):
    if hasattr(x, "dtype"):
        return _type_repr_dtype(x)
    return "Any"


@_type_repr.register(sparse.spmatrix)
@_type_repr.register(np.ndarray)
def _type_repr_dtype(x):
    return escape(str(x.dtype))


@_type_repr.register(pd.DataFrame)
def _type_repr_pd(x: pd.DataFrame):
    return f"DataFrame {escape(str(x.shape))}"


"""
Default Attribute FUNCTIONS
Adds attributes depending on the type of the data structure.
By default does nothing but if overloaded, it can add attributes to the specific data
type automatically.
"""


@singledispatch
def _add_attrs(x, attrs:dict):
    return attrs


@_add_attrs.register(pd.DataFrame)
def _add_attrs_pd_df(x : pd.DataFrame, attrs):
    for c,t in zip(x.columns,x.dtypes):
        txt = None
        if isinstance(t,pd.CategoricalDtype):
            txt = "Ordered " if t.ordered else "Unordered "
            txt += "Cat.: " + ", ".join(t.categories)
            txt = escape(txt)
        else:
            txt = escape(str(t))
        key = escape("col_" + str(c) + "_dtype")
        if attrs.get(key) is None:
            attrs[key] = txt
    return attrs


"""
Inline repr FUNCTIONS
Expected to give short description of the data structure.
"""


@singledispatch
def _inline_format(x, max_width):
    return _maybe_truncate(escape(str(x)), max_width=max_width)


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
    """Summarizes x, gives the html content

    Args:
        name (str): Name to show
        x (Union[pd.DataFrame, np.ndarray, sparse.spmatrix]): The data structure to summarize
        attrs (_type_, optional): Additional information to be
            shown under the attributes button. Defaults to None.

    Returns:
        str: Html repr of x as str
    """

    dims_str = _dim_repr(x)
    name = escape(str(name))
    dtype = _type_repr(x)
    attrs = attrs if attrs else {}
    attrs = _add_attrs(x,attrs)

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())
    disabled = "" if attrs else "disabled"

    preview = _inline_format(x, 35)
    data_repr = _html_format(x)
    attrs_ul = _summarize_attrs(attrs)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    return (
        f"<div class='ad-var-name'><span>{name}</span></div>"
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


def _create_sections_from_conf(x, sections_conf):
    sections = []
    for k, v in sections_conf.items():
        # check if the configuration is filled or not
        if v:
            # TODO: maybe use enums?
            if v["section_type"] == "single":
                sections.append(
                    _collapsible_section(
                        name=k,
                        details=_summarize_item_html(
                            name="", x=getattr(x, k), attrs=v.get("attrs", {})
                        ),
                        **v["args"],
                    )
                )
            elif v["section_type"] == "mapping":
                sections.append(
                    _collapsible_section(
                        name=k,
                        details=_summarize_mapping_html(x=getattr(x, k)),
                        **v["args"],
                    )
                )
            else:
                raise NotImplementedError()

    return sections


def _create_anndata_display_conf(ad_obj: "anndata.AnnData"):
    """Factory to create the configuration of AnnData repr
    the options are hard-coded here (e.g., max_items_collapse options).

    Args:
        ad_obj (anndata.AnnData): data represent

    Returns:
        dict: dict containing the configuration
    """
    max_items_collapse = {
        "X": 1,
        "obs": 1,
        "var": 1,
        "obsm": 1,
        "varm": 1,
        "layers": 1,
        "varp": 1,
        "obsp": 1,
        "uns": 1,
    }

    sections_conf = {
        # sections consisting of single items like matrices
        "X": {
            "section_type": "single",
            "attrs": {
                "Is backed": "Nowhere"
                if not ad_obj.isbacked
                else escape(str(ad_obj.file.filename))
            },
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["X"] < ad_obj.n_obs,
                ),
                "n_items": ad_obj.n_obs,
            },
        }
        if ad_obj.X is not None
        else {},
        "obs": {
            "section_type": "single",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["obs"] < ad_obj.n_obs,
                ),
                "n_items": ad_obj.n_obs,
            },
        }
        if ad_obj.obs is not None
        else {},
        "var": {
            "section_type": "single",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["var"] < ad_obj.n_vars,
                ),
                "n_items": ad_obj.n_vars,
            },
        }
        if ad_obj.var is not None
        else {},
        # dict like sections
        "obsm": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["obsm"] < len(ad_obj.obsm),
                ),
                "n_items": len(ad_obj.obsm),
            }
            if ad_obj.obsm is not None
            else {},
        },
        "varm": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["varm"] < len(ad_obj.varm),
                ),
                "n_items": len(ad_obj.varm),
            }
            if ad_obj.varm is not None
            else {},
        },
        "obsp": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["obsp"] < len(ad_obj.obsp),
                ),
                "n_items": len(ad_obj.obsp),
            }
            if ad_obj.obsp is not None
            else {},
        },
        "varp": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["varp"] < len(ad_obj.varp),
                ),
                "n_items": len(ad_obj.varp),
            }
            if ad_obj.varp is not None
            else {},
        },
        "layers": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["layers"] < len(ad_obj.layers),
                ),
                "n_items": len(ad_obj.layers),
            }
            if ad_obj.layers is not None
            else {},
        },
        "uns": {
            "section_type": "mapping",
            "args": {
                "inline_details": "",
                "enabled": True,
                "collapsed": _get_boolean_with_default(
                    "display_expand_single_item_section",
                    max_items_collapse["uns"] < len(ad_obj.uns),
                ),
                "n_items": len(ad_obj.uns),
            }
            if ad_obj.uns is not None
            else {},
        },
    }
    return sections_conf


def _create_anndata_repr(ad_obj: "anndata.AnnData"):
    obj_type = f"anndata.{type(ad_obj).__name__}"

    header_components = [f"<div class='ad-obj-type'>{escape(obj_type)}</div>"]

    sections_conf = _create_anndata_display_conf(ad_obj)
    sections = _create_sections_from_conf(ad_obj, sections_conf)
    return _obj_repr(ad_obj, header_components, sections)
