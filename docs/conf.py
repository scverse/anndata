import os
import sys
import time
import inspect
from pathlib import Path, PurePosixPath
from typing import Optional, Union, Mapping
import logging
from datetime import datetime

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / '..'))

import anndata


logger = logging.getLogger(__name__)

for generated in HERE.glob('anndata.*.rst'):
    generated.unlink()


# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    # 'plot_generator',
    # 'plot_directive',
    'sphinx_autodoc_typehints',  # needs to be after napoleon
    # 'ipython_directive',
    # 'ipython_console_highlighting',
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_mock_imports = ['_tkinter']  # why this?
autodoc_member_order = 'bysource'
#autodoc_default_flags = ['members']
napoleon_use_rtype = False
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]

intersphinx_mapping = dict(
    python=('https://docs.python.org/3', None),
    h5py=('http://docs.h5py.org/en/latest/', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
    pandas=('http://pandas.pydata.org/pandas-docs/stable/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
)

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'anndata'
author = 'Alex Wolf, Philipp Angerer, Sergei Rybakov'
copyright = f'{datetime.now():%Y}, {author}.'

version = anndata.__version__.replace('.dirty', '')
release = version
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False

# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 2,
}
html_context = dict(
    display_github=True,      # Integrate GitHub
    github_user='theislab',   # Username
    github_repo='anndata',    # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',    # Path in the checkout to the docs root
    css_files=[               # Additional CSS
        '_static/css/custom.css',
    ],
)
if 'READTHEDOCS' in os.environ:
    # For some reason, RTD doesn’t insert their stuff anymore once we add custom CSS files.
    html_context['css_files'].insert(0, 'https://media.readthedocs.org/css/sphinx_rtd_theme.css')

html_static_path = ['_static']

# Output file base name for HTML help builder.
htmlhelp_basename = 'anndatadoc'


# -- Options for other output ---------------------------------------

man_pages = [
    (master_doc, 'anndata', 'anndata Documentation',
     [author], 1)
]
texinfo_documents = [
    (master_doc, 'anndata', 'anndata Documentation',
     author, 'anndata', 'One line description of project.',
     'Miscellaneous'),
]


# -- GitHub links ---------------------------------------------------------


def get_obj_module(qualname):
    """Get a module/class/attribute and its original module by qualname"""
    modname = qualname
    classname = None
    attrname = None
    while modname not in sys.modules:
        attrname = classname
        modname, classname = modname.rsplit('.', 1)

    # retrieve object and find original module name
    if classname:
        cls = getattr(sys.modules[modname], classname)
        modname = cls.__module__
        obj = getattr(cls, attrname) if attrname else cls
    else:
        obj = None

    return obj, sys.modules[modname]


def get_linenos(obj):
    """Get an object’s line numbers"""
    try:
        lines, start = inspect.getsourcelines(obj)
    except TypeError:
        return None, None
    else:
        return start, start + len(lines) - 1


project_dir = Path(__file__).parent.parent  # project/docs/conf.py/../.. → project/
github_url1 = 'https://github.com/{github_user}/{github_repo}/tree/{github_version}'.format_map(html_context)
github_url2 = 'https://github.com/theislab/anndata/tree/master'


def modurl(qualname):
    """Get the full GitHub URL for some object’s qualname."""
    obj, module = get_obj_module(qualname)
    github_url = github_url1
    try:
        path = PurePosixPath(Path(module.__file__).resolve().relative_to(project_dir))
    except ValueError:
        # trying to document something from another package
        github_url = github_url2
        path = '/'.join(module.__file__.split('/')[-2:])
    start, end = get_linenos(obj)
    fragment = '#L{}-L{}'.format(start, end) if start and end else ''
    return '{}/{}{}'.format(github_url, path, fragment)


# html_context doesn’t apply to autosummary templates ☹
# and there’s no way to insert filters into those templates
# so we have to modify the default filters
from jinja2.defaults import DEFAULT_FILTERS

DEFAULT_FILTERS['modurl'] = modurl


# -- Override some classnames in autodoc --------------------------------------------
# This makes sure that automatically documented links actually
# end up being links instead of pointing nowhere.


import sphinx_autodoc_typehints

qualname_overrides = {
    'anndata.base.AnnData': 'anndata.AnnData',
    'pandas.core.frame.DataFrame': 'pandas.DataFrame',
    'scipy.sparse.base.spmatrix': 'scipy.sparse.spmatrix',
    'scipy.sparse.csr.csr_matrix': 'scipy.sparse.csr_matrix',
    'scipy.sparse.csc.csc_matrix': 'scipy.sparse.csc_matrix',
}

fa_orig = sphinx_autodoc_typehints.format_annotation
def format_annotation(annotation):
    # display `Union[A, B]` as `A, B`
    if getattr(annotation, '__origin__', None) is Union or hasattr(annotation, '__union_params__'):
        params = getattr(annotation, '__union_params__', None) or getattr(annotation, '__args__', None)
        # never use the `Optional` keyword in the displayed docs, instead, use the more verbose `, None`
        # as is the convention in the other large numerical packages
        # if len(params or []) == 2 and getattr(params[1], '__qualname__', None) == 'NoneType':
        #     return fa_orig(annotation)  # Optional[...]
        return ', '.join(map(format_annotation, params))
    # do not show the arguments of Mapping
    if getattr(annotation, '__origin__', None) is Mapping:
        return ':class:`~typing.Mapping`'
    if inspect.isclass(annotation):
        full_name = '{}.{}'.format(annotation.__module__, annotation.__qualname__)
        override = qualname_overrides.get(full_name)
        if override is not None:
            return f':py:class:`~{override}`'
    return fa_orig(annotation)
sphinx_autodoc_typehints.format_annotation = format_annotation


# -- Change default role --------------------------------------------

from docutils.parsers.rst import roles

roles.DEFAULT_INTERPRETED_ROLE = 'literal'


# -- Prettier Param docs --------------------------------------------


from typing import Dict, List, Tuple

from docutils import nodes
from sphinx import addnodes
from sphinx.domains.python import PyTypedField, PyObject
from sphinx.environment import BuildEnvironment


class PrettyTypedField(PyTypedField):
    list_type = nodes.definition_list
    
    def make_field(
        self,
        types: Dict[str, List[nodes.Node]],
        domain: str,
        items: Tuple[str, List[nodes.inline]],
        env: BuildEnvironment = None
    ) -> nodes.field:
        def makerefs(rolename, name, node):
            return self.make_xrefs(rolename, domain, name, node, env=env)
        
        def handle_item(fieldarg: str, content: List[nodes.inline]) -> nodes.definition_list_item:
            head = nodes.term()
            head += makerefs(self.rolename, fieldarg, addnodes.literal_strong)
            fieldtype = types.pop(fieldarg, None)
            if fieldtype is not None:
                head += nodes.Text(' : ')
                if len(fieldtype) == 1 and isinstance(fieldtype[0], nodes.Text):
                    text_node, = fieldtype  # type: nodes.Text
                    head += makerefs(self.typerolename, text_node.astext(), addnodes.literal_emphasis)
                else:
                    head += fieldtype

            body_content = nodes.paragraph('', '', *content)
            body = nodes.definition('', body_content)

            return nodes.definition_list_item('', head, body)

        fieldname = nodes.field_name('', self.label)
        if len(items) == 1 and self.can_collapse:
            fieldarg, content = items[0]
            bodynode = handle_item(fieldarg, content)
        else:
            bodynode = self.list_type()
            for fieldarg, content in items:
                bodynode += handle_item(fieldarg, content)
        fieldbody = nodes.field_body('', bodynode)
        return nodes.field('', fieldname, fieldbody)


# replace matching field types with ours
PyObject.doc_field_types = [
    PrettyTypedField(
        ft.name,
        names=ft.names,
        typenames=ft.typenames,
        label=ft.label,
        rolename=ft.rolename,
        typerolename=ft.typerolename,
        can_collapse=ft.can_collapse,
    ) if isinstance(ft, PyTypedField) else ft
    for ft in PyObject.doc_field_types
]

