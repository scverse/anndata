import sys

if sys.version_info < (3, 6):
    sys.exit("anndata requires Python >= 3.6")
from pathlib import Path

from setuptools import setup, find_namespace_packages

try:
    import pytoml
except ImportError:
    sys.exit("Please use `pip install .` or install pytoml first.")

proj = pytoml.loads(Path("pyproject.toml").read_text())
metadata = proj["tool"]["anndata"]

setup(
    name="anndata",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    description="Annotated Data.",
    long_description=Path("README.rst").read_text("utf-8"),
    url="http://github.com/theislab/anndata",
    author=metadata["author"],
    author_email=metadata["author-email"],
    license="BSD-3-Clause",
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    extras_require=dict(
        dev=["setuptools_scm", "pytoml", "black>=20.8b1"],
        doc=[
            # Sphinx 2 has nicer looking sections
            "sphinx>=2.0.1",
            "sphinx-rtd-theme",
            "sphinx-autodoc-typehints>=1.11.0",
            "sphinx_issues",
            "scanpydoc>=0.5",
            'typing_extensions; python_version < "3.8"',
        ],
        test=[
            "loompy>=3.0.5",
            "pytest>=4.6",
            "pytest-cov>=2.10",
            "codacy-coverage",
            "docutils",  # for rst2html.py
            "zarr",
            "matplotlib",
            "sklearn",
            "xlrd",
            "joblib",
            "boltons",
            "scanpy",
        ],
    ),
    python_requires=">=3.6",
    packages=find_namespace_packages(include=["anndata", "anndata.*"]),
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Environment :: Console",
        "Framework :: Jupyter",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
