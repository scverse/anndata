from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path


def pytest_ignore_collect(collection_path: Path) -> bool:
    return collection_path.name == "hv.py"
