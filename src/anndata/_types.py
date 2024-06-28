"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import Union

from anndata.compat import H5Array, H5Group, ZarrArray, ZarrGroup

__all__ = [
    "ArrayStorageType",
    "GroupStorageType",
    "StorageType",
]

ArrayStorageType = Union[ZarrArray, H5Array]
GroupStorageType = Union[ZarrGroup, H5Group]
StorageType = Union[ArrayStorageType, GroupStorageType]
