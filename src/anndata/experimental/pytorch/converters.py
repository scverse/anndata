"""Helper converters for AnnLoader batches.

This module provides convenience converters that can be passed to the
``batch_converter`` parameter of :pyclass:`~anndata.experimental.pytorch.AnnLoader`.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd
import torch

try:  # keep mypy happy when torch not present for docs build
    from torch import Tensor
except ImportError:  # pragma: no cover
    Tensor = Any  # type: ignore

__all__ = ["to_tensor_dict"]


def _to_tensor(arr) -> Tensor | Any:
    """Best-effort conversion of *arr* to ``torch.Tensor``.

    Falls back to returning *arr* unchanged if torch or numpy is not available.
    """
    if isinstance(arr, torch.Tensor):
        return arr
    try:
        import numpy as np
        from scipy.sparse import issparse

        if issparse(arr):
            arr = arr.toarray()
        if isinstance(arr, (np.ndarray, list)):
            return torch.tensor(arr)
    except ImportError:  # pragma: no cover
        pass
    return arr


def to_tensor_dict(batch: Any) -> dict[str, Any]:
    """Convert an AnnLoader batch to a plain ``dict`` of tensors/arrays.

    * ``X``   → ``"x"``
    * each column in ``obs`` becomes a key in the output dict
    * if *batch* is already a mapping it is returned as a *shallow copy*.
    """
    # If user already returns dict-like we preserve it
    if isinstance(batch, Mapping):
        return dict(batch)

    out: dict[str, Any] = {}

    # AnnCollectionView has .X and .obs attributes
    if hasattr(batch, "X"):
        out["x"] = _to_tensor(batch.X)

    if hasattr(batch, "obs") and batch.obs is not None:
        obs_data = batch.obs
        # Handle pandas DataFrame
        if isinstance(obs_data, pd.DataFrame):
            for col in obs_data.columns:
                # ensure unique keys – users can post-process if needed
                out[col] = _to_tensor(obs_data[col].to_numpy())
        # Handle AnnCollection MapObsView (can be converted to dict directly)
        elif hasattr(obs_data, "to_dict"):
            obs_dict = obs_data.to_dict()
            for key, value in obs_dict.items():
                out[key] = _to_tensor(value)
        # Handle generic dict-like objects
        elif hasattr(obs_data, "keys") and callable(obs_data.keys):
            try:
                obs_dict = dict(obs_data)
                for key, value in obs_dict.items():
                    out[key] = _to_tensor(value)
            except (TypeError, AttributeError, ValueError):
                pass  # Skip if conversion fails

    return out
