from __future__ import annotations

from typing import TYPE_CHECKING, overload

import pandas as pd
from holoviews.core.dimension import Dimension

from . import AdAcc, AdPath, GraphVecAcc, LayerVecAcc, MetaVecAcc, MultiVecAcc

if TYPE_CHECKING:
    from typing import Any, Literal, Self

    from . import VecAcc


__all__ = ["A", "AdDim"]


def mk_label(p: AdPath, /) -> str | None:
    match p.acc:
        case MultiVecAcc():
            return f"{p.acc.k} {p.idx}"
        case GraphVecAcc():
            return next((f"{p.acc.k} {i}" for i in p.idx if isinstance(i, str)), None)
        case LayerVecAcc():
            return next(
                (
                    f"{p.acc.k} {i}" if p.acc.k else i
                    for i in p.idx
                    if isinstance(i, str)
                ),
                None,
            )
        case MetaVecAcc():
            return f"{p.acc.ax} index" if p.idx is pd.Index else p.idx
        case _:  # pragma: no cover
            msg = f"Unsupported vector accessor {p.acc!r}"
            raise AssertionError(msg)


class AdDim[I](AdPath[I], Dimension):
    def __init__(
        self,
        acc: VecAcc[Self, I],
        idx: I,
        /,
        *,
        label: str | None = None,
        **params: object,
    ) -> None:
        label_kw = {} if label is None else dict(label=label)
        AdPath.__init__(self, acc, idx)  # type: ignore
        Dimension.__init__(self, repr(self), label=mk_label(self), **params, **label_kw)

    @overload
    @classmethod
    def from_dimension(
        cls, dim: Dimension, *, strict: Literal[True] = True
    ) -> Self: ...
    @overload
    @classmethod
    def from_dimension(
        cls, dim: Dimension, *, strict: Literal[False]
    ) -> Self | None: ...

    @classmethod
    def from_dimension(cls, dim: Dimension, *, strict: bool = True) -> Self | None:
        """Create accessor from another dimension."""
        if TYPE_CHECKING:
            assert isinstance(dim.name, str)

        if isinstance(dim, cls):
            return dim
        if (rv := A.resolve(dim.name, strict=strict)) is None:
            return None
        if dim.name != dim.label:
            rv.label = dim.label
        return rv

    def clone(
        self,
        spec: str | tuple[str, str] | None = None,
        **overrides: Any,
    ) -> Self:
        """Clones the Dimension with new parameters.

        Derive a new Dimension that inherits existing parameters
        except for the supplied, explicit overrides

        Parameters
        ----------
        spec
            Dimension tuple specification
        func
            Function to resolve the dimension values
            given the AnnData object.
        axes
            The axes represented by the Dimension
        overrides
            Dimension parameter overrides

        Returns
        -------
        Cloned Dimension object

        """
        settings = dict(self.param.values(), **overrides)
        func = settings.pop("func", self._func)
        axes = settings.pop("axes", self.axes)

        if spec is None:
            spec = (self.name, overrides.get("label", self.label))
        if "label" in overrides and isinstance(spec, str):
            spec = (spec, overrides["label"])
        elif "label" in overrides and isinstance(spec, tuple):
            if overrides["label"] != spec[1]:
                self.param.warning(
                    f"Using label as supplied by keyword ({overrides['label']!r}), "
                    f"ignoring tuple value {spec[1]!r}"
                )
            spec = (spec[0], overrides["label"])
        return self.__class__(
            spec,
            func,
            axes,
            **{k: v for k, v in settings.items() if k not in ["name", "label"]},
        )

    def __hash__(self) -> int:
        return hash(self._repr)

    def __eq__(self, dim: object) -> bool:
        # shortcut if label, number, or so matches
        if super().__eq__(dim):
            return True
        # try to resolve
        if isinstance(dim, str) and (dim := A.resolve(dim, strict=False)) is None:
            return False
        # if dim is a non-matching dimension (e.g. from a string), convert
        if isinstance(dim, Dimension):
            if not isinstance(dim, AdPath):
                if dim.name == self.name:
                    return True
                if (dim := type(self).from_dimension(dim, strict=False)) is None:
                    return False
            # dim is an AdPath, check equality
            return hash(self) == hash(dim)
        # some unknown type
        return False


A = AdAcc(path_class=AdDim)
