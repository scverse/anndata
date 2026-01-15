from __future__ import annotations

from typing import TYPE_CHECKING, overload

import pandas as pd
from holoviews.core.dimension import Dimension

from . import AdAcc, AdRef, GraphVecAcc, LayerVecAcc, MetaVecAcc, MultiVecAcc

if TYPE_CHECKING:
    from typing import Any, Literal, Self

    from . import VecAcc


__all__ = ["A", "AdDim"]


def mk_label[I](p: AdRef[I], /) -> str | None:
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


class AdDim[I](AdRef[I], Dimension):
    def __init__(self, acc: VecAcc[Self, I], idx: I, /, **params: object) -> None:
        AdRef.__init__(self, acc, idx)  # type: ignore
        spec = params.pop("spec", repr(self))
        label = params.pop("label", mk_label(self))
        Dimension.__init__(self, spec, label=label, **params)

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
        overrides
            Dimension parameter overrides

        Returns
        -------
        Cloned Dimension object

        """
        settings = dict(self.param.values(), **overrides)
        acc = settings.pop("acc", self.acc)
        idx = settings.pop("idx", self.idx)

        match spec, ("label" in overrides):
            case None | str(), _:
                spec = (spec or self.name, overrides.get("label", self.label))
            case (name, label), True:
                if overrides["label"] != label:
                    self.param.warning(
                        f"Using label as supplied by keyword ({overrides['label']!r}), "
                        f"ignoring tuple value {label!r}"
                    )
                spec = (name, overrides["label"])

        return type(self)(
            acc,
            idx,
            **{k: v for k, v in settings.items() if k not in ["name", "label"]},
        )

    def __hash__(self) -> int:
        return hash((type(self), repr(self)))

    def __eq__(self, dim: object) -> bool:
        # shortcut if label, number, or so matches
        if super().__eq__(dim):
            return True
        # try to resolve
        if isinstance(dim, str) and (dim := A.resolve(dim, strict=False)) is None:
            return False
        # if dim is a non-matching dimension (e.g. from a string), convert
        if isinstance(dim, Dimension):
            if not isinstance(dim, AdRef):
                if dim.name == self.name:
                    return True
                if (dim := type(self).from_dimension(dim, strict=False)) is None:
                    return False
            # dim is an AdRef, check equality
            return hash(self) == hash(dim)
        # some unknown type
        return False


A = AdAcc(path_class=AdDim)
