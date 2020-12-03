from pathlib import Path

here = Path(__file__).parent

try:
    from setuptools_scm import get_version
    import pytoml

    proj = pytoml.loads((here.parent / "pyproject.toml").read_text())
    metadata = proj["tool"]["anndata"]

    __version__ = get_version(root="..", relative_to=__file__)
    __author__ = metadata["author"]
    __email__ = metadata["author-email"]
except (ImportError, LookupError, FileNotFoundError):
    try:
        from importlib.metadata import metadata
    except ImportError:
        from importlib_metadata import metadata

    meta = metadata(here.name)
    __version__ = meta["Version"]
    __author__ = meta["Author"]
    __email__ = meta["Author-email"]
