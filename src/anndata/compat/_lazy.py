"""Stand-ins for classes from optional dependencies that don’t import them.

Importing dask, xarray or awkward costs more than importing anndata itself,
so :mod:`anndata.compat` exposes a :class:`LazyType` for each of their classes instead.
A lazy type behaves like the class it stands for in

- :func:`isinstance` and :func:`issubclass` (and therefore in ``match`` class patterns),
- :func:`functools.singledispatch`, which treats it as a virtual base class, and
- calls and attribute access (e.g. ``XDataset.from_dataframe(...)``).

The first two never import anything: an object can only be an instance of a class
whose module has been imported, so until the *guard* module is in :data:`sys.modules`,
the answer is simply “no”.
Calls and attribute access do import the module, and fail if it isn’t installed.

Mappings keyed by type (such as the IO registry) need to look up
:func:`canonical_type` of a type, since a lazy type and its real class are different objects.
"""

from __future__ import annotations

import functools
import sys
from importlib import import_module
from types import UnionType
from typing import TYPE_CHECKING, cast, get_args, get_type_hints

if TYPE_CHECKING:
    from collections.abc import Callable
    from functools import _SingleDispatchCallable as _SingleDispatch
    from typing import Any


__all__ = ["LazyType", "canonical_type", "lazy_type", "resolve", "singledispatch"]


_PROXY_FOR_REAL: dict[type, LazyType] = {}
_UNRESOLVED: dict[LazyType, None] = {}  # ordered set


class LazyType(type):
    """Metaclass of a stand-in for a class that is only imported once needed."""

    _lazy_guard: str
    _lazy_load: Callable[[], type]
    _lazy_real: type | None

    def _lazy_resolve(cls, *, force: bool = False) -> type | None:
        """Return the real class, or `None` if its module hasn’t been imported (and `force` is false)."""
        if (real := cls.__dict__["_lazy_real"]) is not None:
            return real
        if not force and cls._lazy_guard not in sys.modules:
            return None
        try:
            real = cls.__dict__["_lazy_load"]()
        except (ImportError, AttributeError):
            # e.g. the guard module is still being initialized
            if force:
                raise
            return None
        type.__setattr__(cls, "_lazy_real", real)
        _PROXY_FOR_REAL[real] = cls
        _UNRESOLVED.pop(cls, None)
        return real

    def _lazy_force(cls) -> type:
        """Return the real class, importing its module if necessary."""
        real = cls._lazy_resolve(force=True)
        assert real is not None
        return real

    def __instancecheck__(cls, obj: object) -> bool:
        real = cls._lazy_resolve()
        return real is not None and isinstance(obj, real)

    def __subclasscheck__(cls, sub: type) -> bool:
        if isinstance(sub, LazyType):
            if cls in sub.__mro__:
                return True
            if (sub_real := sub._lazy_resolve()) is None:
                return False
            sub = sub_real
        real = cls._lazy_resolve()
        return real is not None and issubclass(sub, real)

    def __call__(cls, *args: object, **kwargs: object) -> object:
        return cls._lazy_force()(*args, **kwargs)

    def __getattr__(cls, name: str) -> object:
        # Dunder lookups come from introspection (e.g. `singledispatch` checking for
        # `__abstractmethods__`) and must not trigger an import.
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        return getattr(cls._lazy_force(), name)

    def __repr__(cls) -> str:
        return f"<lazy class '{cls.__module__}.{cls.__qualname__}'>"


def lazy_type(
    module: str,
    qualname: str,
    *,
    guard: str | None = None,
    load: Callable[[], type] | None = None,
    bases: tuple[LazyType, ...] = (),
) -> LazyType:
    """Create a stand-in for `module.qualname`.

    Parameters
    ----------
    module
        Module the class is importable from.
    qualname
        Name of the class within `module`.
    guard
        Module whose presence in :data:`sys.modules` means the class may have instances.
        Defaults to `module`.
    load
        Returns the real class. Defaults to importing `module.qualname`.
        Use this for subclasses of optional classes that anndata defines itself.
    bases
        Lazy types the real class subclasses.
        Needed so :func:`functools.singledispatch` knows which registration is more specific.
    """
    if load is None:

        def load() -> type:
            obj: Any = import_module(module)
            for part in qualname.split("."):
                obj = getattr(obj, part)
            return obj

    cls = LazyType(
        qualname.rpartition(".")[2],
        bases,
        dict(
            __module__=module,
            __qualname__=qualname,
            _lazy_guard=guard or module,
            _lazy_load=staticmethod(load),
            _lazy_real=None,
        ),
    )
    _UNRESOLVED[cls] = None
    return cls


def singledispatch[R](func: Callable[..., R]) -> _SingleDispatch[R]:
    """:func:`functools.singledispatch` that accepts :class:`LazyType`\\ s.

    A lazy type is registered as its real class once that class’s module has been imported,
    so dispatch then resolves exactly as if the class had been registered directly.
    (Registering the lazy type itself would make it a *virtual* base class,
    which loses against real ABC bases, e.g. ``xarray.Dataset`` is a real :class:`~collections.abc.Mapping`.)
    """
    dispatcher = functools.singledispatch(func)
    pending: dict[LazyType, Callable[..., Any]] = {}
    n_modules = -1

    def resolve_pending() -> None:
        nonlocal n_modules
        # cheap check: a lazy type can only become resolvable after an import
        if not pending or n_modules == len(sys.modules):
            return
        n_modules = len(sys.modules)
        for proxy in list(pending):
            if (real := proxy._lazy_resolve()) is not None:
                dispatcher.register(real, pending.pop(proxy))

    def register(
        cls: type | UnionType | Callable[..., Any],
        impl: Callable[..., Any] | None = None,
    ) -> Callable[..., Any]:
        if impl is None:
            if isinstance(cls, type | UnionType):
                return lambda f: register(cls, f)
            # `@f.register` with a type annotation on the first parameter
            impl = cls
            cls = next(iter(get_type_hints(impl).values()))
        if isinstance(cls, UnionType):
            for c in get_args(cls):
                register(c, impl)
            return impl
        if isinstance(cls, LazyType):
            if (real := cls._lazy_resolve()) is None:
                nonlocal n_modules
                pending[cls] = impl
                n_modules = -1
                return impl
            cls = real
        dispatcher.register(cast("type", cls), impl)
        return impl

    def dispatch(cls: type) -> Callable[..., R]:
        if isinstance(cls, LazyType):
            cls = cls._lazy_force()
        resolve_pending()
        return dispatcher.dispatch(cls)

    @functools.wraps(func)
    def wrapper(*args: object, **kwargs: object) -> R:
        if not args:
            msg = f"{funcname} requires at least 1 positional argument"
            raise TypeError(msg)
        return dispatch(args[0].__class__)(*args, **kwargs)

    funcname = getattr(func, "__name__", "singledispatch function")
    wrapper.register = register  # type: ignore[attr-defined]
    wrapper.dispatch = dispatch  # type: ignore[attr-defined]
    wrapper.registry = dispatcher.registry  # type: ignore[attr-defined]
    wrapper._clear_cache = dispatcher._clear_cache  # type: ignore[attr-defined]
    return cast("_SingleDispatch[R]", wrapper)


def resolve(typ: type) -> type:
    """Return the real class of a :class:`LazyType` (importing it), or `typ` itself."""
    return typ._lazy_force() if isinstance(typ, LazyType) else typ


def canonical_type[T](typ: T) -> T | LazyType:
    """Map a real class to its :class:`LazyType` stand-in, if it has one."""
    if (proxy := _PROXY_FOR_REAL.get(typ)) is not None:  # type: ignore[call-overload]
        return proxy
    if not _UNRESOLVED or isinstance(typ, LazyType) or not isinstance(typ, type):
        return typ
    # the class might be from a module that was imported after we last looked
    for proxy in list(_UNRESOLVED):
        proxy._lazy_resolve()
    return _PROXY_FOR_REAL.get(typ, typ)
