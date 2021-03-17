from typing import Iterable
import math
import re
import numpy as np

def make_single_column_table(heading: str, entries: Iterable[str]) -> str:
    table = f"""
        <table>
            <thead>
                <tr><th> {heading} </th></tr>
            </thead>
            <tbody>
            """
    for entry in entries:
        table += f"""
                <tr>
                    <td>{entry}</td>
                </tr>
                """
    table +=  """
            </tbody>
        </table>
    """  
    return table

def make_single_row_table(entries: Iterable[str]) -> str:
    table = f"""
        <table>
            <tbody>
                <tr>
            """
    for entry in entries:
        table += f"""<td style="vertical-align:top">{entry}</td>"""
    table +=  """
                </tr>
            </tbody>
        </table>
    """
    return table

# much of the svg code below is copied or adapted
# from the dask project 
# https://github.com/dask/dask/blob/main/dask/array/svg.py


text_style = 'font-size="1.0rem" font-weight="100" text-anchor="middle"'


def svg_anndata(n_obs, n_var, size=200, sizes=None):
    shape = (n_var, n_obs)
    sizes = sizes or draw_sizes(shape, size=size)
    lines, (min_x, max_x, min_y, max_y) = svg_box(sizes[0], sizes[1], size=size)

    header = (
        '<svg width="%d" height="%d" style="stroke:rgb(0,0,0);stroke-width:1" >\n'
        % (max_x + 50, max_y + 50)
    )
    footer = "\n</svg>"

    #if shape[0] >= 100:
    #    rotate = -90
    #else:
    rotate = 0

    text = [
        "",
        "  <!-- Text -->",
        '  <text x="%f" y="%f" %s >%d</text>'
        % (max_x / 2, max_y + 20, text_style, shape[0]),
        '  <text x="%f" y="%f" %s transform="rotate(%d,%f,%f)">%d</text>'
        % (max_x + 20, max_y / 2, text_style, rotate, max_x + 20, max_y / 2, shape[1]),
    ]

    return header + "\n".join(lines + text) + footer



def svg_lines(x1, y1, x2, y2):
    """Convert points into lines of text for an SVG plot
    Examples
    --------
    >>> svg_lines([0, 1], [0, 0], [10, 11], [1, 1])  # doctest: +NORMALIZE_WHITESPACE
    ['  <line x1="0" y1="0" x2="10" y2="1" style="stroke-width:2" />',
     '  <line x1="1" y1="0" x2="11" y2="1" style="stroke-width:2" />']
    """
    n = len(x1)

    lines = [
        '  <line x1="%d" y1="%d" x2="%d" y2="%d" />' % (x1[i], y1[i], x2[i], y2[i])
        for i in range(n)
    ]

    lines[0] = lines[0].replace(" /", ' style="stroke-width:2" /')
    lines[-1] = lines[-1].replace(" /", ' style="stroke-width:2" /')
    return lines


def svg_box(w, h, size=200):
    """Create lines of SVG text that show a grid
    Parameters
    ----------
    x: numpy.ndarray
    y: numpy.ndarray
    offset: tuple
        translational displacement of the grid in SVG coordinates
    """

    x = np.array((0,w))
    y = np.array((0,h))
    
    # Horizontal lines
    x1 = np.zeros_like(y)
    y1 = y
    x2 = np.full_like(y, x[-1])
    y2 = y

    h_lines = ["", "  <!-- Horizontal lines -->"] + svg_lines(x1, y1, x2, y2)

    # Vertical lines
    x1 = x 
    y1 = np.zeros_like(x) 
    x2 = x 
    y2 = np.full_like(x, y[-1])

    v_lines = ["", "  <!-- Vertical lines -->"] + svg_lines(x1, y1, x2, y2)

    color = "ECB172" # if len(x) < max_n and len(y) < max_n else "8B4903"
    corners = f"{x1[0]},{y1[0]} {x1[-1]},{y1[-1]} {x2[-1]},{y2[-1]} {x2[0]},{y2[0]}"
    rect = [
        "",
        "  <!-- Colored Rectangle -->",
        f'  <polygon points="{corners}" style="fill:#{color}A0;stroke-width:0"/>',
    ]

    return h_lines + v_lines + rect, (0, w, 0, h)


def draw_sizes(shape, size=200):
    """ Get size in pixels for all dimensions """
    mx = max(shape)
    ratios = [mx / max(0.1, d) for d in shape]
    ratios = [ratio_response(r) for r in ratios]
    return tuple(size / r for r in ratios)


def ratio_response(x):
    """How we display actual size ratios
    Common ratios in sizes span several orders of magnitude,
    which is hard for us to perceive.
    We keep ratios in the 1-3 range accurate, and then apply a logarithm to
    values up until about 100 or so, at which point we stop scaling.
    """
    if x < math.e:
        return x
    elif x <= 100:
        return math.log(x + 12.4)  # f(e) == e
    else:
        return math.log(100 + 12.4)