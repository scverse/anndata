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
    >>> class MyTransform(Transform):
    ...     def __init__(self, param=1.0):
    ...         self.param = param
    ...
    ...     def __call__(self, x):
    ...         return x * self.param
    ...
    ...     def __repr__(self):
    ...         return f"MyTransform(param={self.param})"
    """

    @abstractmethod
    def __call__(self, x: torch.Tensor) -> torch.Tensor:
        """Apply the transform to the input tensor."""

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
    >>> class LogTransform(Transform):
    ...     def __call__(self, x):
    ...         return torch.log1p(x)
    >>>
    >>> class NormalizeTransform(Transform):
    ...     def __init__(self, target_sum=1e4):
    ...         self.target_sum = target_sum
    ...
    ...     def __call__(self, x):
    ...         row_sum = torch.sum(x, dim=-1, keepdim=True) + 1e-8
    ...         return x * (self.target_sum / row_sum)
    >>>
    >>> transform = Compose([
    ...     NormalizeTransform(target_sum=1e4),
    ...     LogTransform(),
    ... ])
    """

    def __init__(self, transforms: Sequence[Transform]):
        if not transforms:
            msg = "At least one transform must be provided"
            raise ValueError(msg)

        self.transforms = list(transforms)

    def __call__(self, x: torch.Tensor) -> torch.Tensor:
        """Apply all transforms in sequence."""
        for transform in self.transforms:
            x = transform(x)
        return x

    def __repr__(self) -> str:
        format_string = f"{self.__class__.__name__}("
        for t in self.transforms:
            format_string += f"\n    {t}"
        format_string += "\n)"
        return format_string
