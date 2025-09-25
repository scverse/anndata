"""Transform base classes for AnnDataset following torchvision patterns.

This module provides the base Transform class and Compose functionality
that can be serialized properly for multiprocessing, avoiding pickling
issues with lambda functions and closures.

Users can inherit from Transform to create their own serializable transforms.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

try:
    import torch

    TORCH_AVAILABLE = True
except ImportError:
    # Mock torch for when it's not available
    class torch:
        class Tensor:
            pass

    TORCH_AVAILABLE = False


class Transform(ABC):
    """Abstract base class for all transforms.

    Inherit from this class to create custom transforms that work
    seamlessly with multiprocessing.

    Examples
    --------
    >>> class MyTransform(Transform):  # doctest: +SKIP
    ...     def __init__(self, param=1.0):  # doctest: +SKIP
    ...         self.param = param  # doctest: +SKIP
    ...
    ...     def __call__(self, data_dict):  # doctest: +SKIP
    ...         data_dict["X"] = data_dict["X"] * self.param  # doctest: +SKIP
    ...         return data_dict  # doctest: +SKIP
    ...
    ...     def __repr__(self):  # doctest: +SKIP
    ...         return f"MyTransform(param={self.param})"  # doctest: +SKIP
    """

    @abstractmethod
    def __call__(self, data_dict: dict) -> dict:
        """Apply the transform to the input data dictionary.

        Parameters
        ----------
        data_dict : dict
            Complete data dictionary with 'X' and metadata keys (obs_*)

        Returns
        -------
        dict
            Transformed data dictionary
        """

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"


class Compose(Transform):
    """Compose several transforms together.

    Similar to torchvision.transforms.Compose, this allows chaining
    multiple transforms while remaining serializable for multiprocessing.

    Parameters
    ----------
    transforms
        List of transform objects to compose

    Examples
    --------
    >>> class LogTransform(Transform):  # doctest: +SKIP
    ...     def __call__(self, data_dict):  # doctest: +SKIP
    ...         data_dict["X"] = torch.log1p(data_dict["X"])  # doctest: +SKIP
    ...         return data_dict  # doctest: +SKIP
    >>>
    >>> class NormalizeTransform(Transform):  # doctest: +SKIP
    ...     def __init__(self, target_sum=1e4):  # doctest: +SKIP
    ...         self.target_sum = target_sum  # doctest: +SKIP
    ...
    ...     def __call__(self, data_dict):  # doctest: +SKIP
    ...         X = data_dict["X"]  # doctest: +SKIP
    ...         row_sum = torch.sum(X, dim=-1, keepdim=True) + 1e-8  # doctest: +SKIP
    ...         data_dict["X"] = X * (self.target_sum / row_sum)  # doctest: +SKIP
    ...         return data_dict  # doctest: +SKIP
    >>>
    >>> transform = Compose([  # doctest: +SKIP
    ...     NormalizeTransform(target_sum=1e4),  # doctest: +SKIP
    ...     LogTransform(),  # doctest: +SKIP
    ... ])  # doctest: +SKIP
    """

    def __init__(self, transforms: Sequence[Transform]):
        if not transforms:
            msg = "At least one transform must be provided"
            raise ValueError(msg)

        self.transforms = list(transforms)

    def __call__(self, data_dict: dict) -> dict:
        """Apply all transforms in sequence."""
        for transform in self.transforms:
            data_dict = transform(data_dict)
        return data_dict

    def __repr__(self) -> str:
        format_string = f"{self.__class__.__name__}("
        for t in self.transforms:
            format_string += f"\n    {t}"
        format_string += "\n)"
        return format_string
