from __future__ import annotations

from importlib.util import find_spec

# public re-exports
from ._annloader import AnnLoader

__all__: list[str] = ["AnnLoader"]

# Only import batch_dict_converter if torch is available
if find_spec("torch"):
    from .converters import to_tensor_dict as batch_dict_converter

    __all__ += ["batch_dict_converter"]
else:
    # Provide a fallback that raises a helpful error
    def batch_dict_converter(*args, **kwargs):
        msg = "batch_dict_converter requires PyTorch. Install with: pip install torch"
        raise ImportError(msg)

    __all__ += ["batch_dict_converter"]
