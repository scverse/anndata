"""Get version from VCS in a dev environment or from package metadata in production.

See <https://github.com/maresb/hatch-vcs-footgun-example>.
"""

from __future__ import annotations

from pathlib import Path

from hatchling.metadata.core import ProjectMetadata
from hatchling.plugin.manager import PluginManager
from hatchling.utils.fs import locate_file

__all__ = ["__version__"]


def _get_version_from_vcs() -> str:
    if (pyproject_toml := locate_file(__file__, "pyproject.toml")) is None:
        msg = "pyproject.toml not found although hatchling is installed"
        raise LookupError(msg)
    root = Path(pyproject_toml).parent
    metadata = ProjectMetadata(root=str(root), plugin_manager=PluginManager())
    # Version can be either statically set in pyproject.toml or computed dynamically:
    return metadata.core.version or metadata.hatch.version.cached


try:
    __version__ = _get_version_from_vcs()
except (ImportError, LookupError):
    import importlib.metadata

    __version__ = importlib.metadata.version("anndata")
