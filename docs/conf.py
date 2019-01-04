import sys
import logging
from pathlib import Path
from datetime import datetime

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
import anndata  # noqa


logger = logging.getLogger(__name__)

for generated in HERE.glob('anndata.*.rst'):
    generated.unlink()


# -- General configuration ------------------------------------------------


needs_sphinx = '1.7'  # autosummary bugfix

# General information
project = 'anndata'
author = ', '.join([
    'Philipp Angerer',
    'Alex Wolf',
    'Sergei Rybakov',
])
copyright = f'{datetime.now():%Y}, {author}.'
version = anndata.__version__.replace('.dirty', '')
release = version

# default settings
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
default_role = 'literal'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'

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
    'scanpydoc',
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False

intersphinx_mapping = dict(
    h5py=('http://docs.h5py.org/en/latest/', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    pandas=('http://pandas.pydata.org/pandas-docs/stable/', None),
    python=('https://docs.python.org/3', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
)


# -- Options for HTML output ----------------------------------------------


html_theme = 'sphinx_rtd_theme'
html_theme_options = dict(
    navigation_depth=2,
)
html_context = dict(
    display_github=True,      # Integrate GitHub
    github_user='theislab',   # Username
    github_repo='anndata',    # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',    # Path in the checkout to the docs root
)
html_static_path = ['_static']


def setup(app):
    app.add_stylesheet('css/custom.css')


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f'{project}doc'
doc_title = f'{project} Documentation'
latex_documents = [
    (master_doc, f'{project}.tex', doc_title, author, 'manual'),
]
man_pages = [
    (master_doc, project, doc_title, [author], 1)
]
texinfo_documents = [
    (master_doc, project, doc_title, author, project, 'One line description of project.', 'Miscellaneous'),
]
