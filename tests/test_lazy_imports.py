"""Optional dependencies are only imported once used (https://github.com/scverse/anndata/issues/1122)."""

from __future__ import annotations

import subprocess
import sys
from abc import ABC
from collections.abc import Mapping
from types import ModuleType
from typing import TYPE_CHECKING

import pytest

from anndata.compat import LazyType, canonical_type, lazy_type, singledispatch

if TYPE_CHECKING:
    from collections.abc import Generator


OPTIONAL_DEPS = ["awkward", "dask", "xarray", "zappy"]


def test_import_does_not_import_optional_deps() -> None:
    code = (
        "import sys, anndata, anndata.experimental; "
        f"print(*sorted(m for m in {OPTIONAL_DEPS!r} if m in sys.modules))"
    )
    loaded = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    ).stdout.split()
    assert loaded == []


@pytest.fixture
def fake_module() -> Generator[str]:
    """Name of a module that is only “imported” once `install_fake_module` is called."""
    name = "_anndata_test_lazy_fake_module"
    yield name
    sys.modules.pop(name, None)


def install_fake_module(name: str) -> tuple[type, type]:
    class Thing(ABC):  # noqa: B024
        kind = "thing"

    class SubThing(Thing):
        pass

    mod = ModuleType(name)
    mod.__dict__.update(Thing=Thing, SubThing=SubThing)
    sys.modules[name] = mod
    return Thing, SubThing


def test_lazy_type_before_import(fake_module: str) -> None:
    Thing = lazy_type(fake_module, "Thing")
    assert isinstance(Thing, LazyType)
    assert not isinstance(1, Thing)
    assert not issubclass(int, Thing)
    assert fake_module not in sys.modules


def test_lazy_type_after_import(fake_module: str) -> None:
    LazyThing = lazy_type(fake_module, "Thing")
    LazySubThing = lazy_type(fake_module, "SubThing", bases=(LazyThing,))
    Thing, SubThing = install_fake_module(fake_module)

    assert isinstance(SubThing(), LazyThing)
    assert isinstance(SubThing(), LazySubThing)
    assert not isinstance(Thing(), LazySubThing)
    assert issubclass(SubThing, LazyThing)
    assert issubclass(LazySubThing, LazyThing)
    assert not issubclass(LazyThing, LazySubThing)
    assert isinstance(SubThing(), LazyThing | int)
    match SubThing():
        case LazyThing():
            pass
        case _:
            pytest.fail("class pattern didn’t match")

    # calls and attributes are forwarded
    assert type(LazySubThing()) is SubThing
    assert LazySubThing.kind == "thing"

    assert canonical_type(Thing) is LazyThing
    assert canonical_type(SubThing) is LazySubThing
    assert canonical_type(int) is int


def test_lazy_type_call_imports(fake_module: str) -> None:
    Thing = lazy_type(fake_module, "Thing")
    with pytest.raises(ModuleNotFoundError):
        Thing()


@pytest.mark.parametrize("register_before_import", [True, False])
def test_singledispatch(fake_module: str, *, register_before_import: bool) -> None:
    @singledispatch
    def f(x: object) -> str:
        return "object"

    LazyThing = lazy_type(fake_module, "Thing")
    LazySubThing = lazy_type(fake_module, "SubThing", bases=(LazyThing,))

    def register() -> None:
        f.register(Mapping, lambda x: "mapping")
        f.register(LazyThing, lambda x: "thing")
        f.register(LazySubThing | int, lambda x: "subthing or int")

    if register_before_import:
        register()
        assert f(1.0) == "object"
        assert f(1) == "subthing or int"
        Thing, SubThing = install_fake_module(fake_module)
    else:
        Thing, SubThing = install_fake_module(fake_module)
        register()

    assert f(Thing()) == "thing"
    assert f(SubThing()) == "subthing or int"
    assert f.dispatch(LazyThing)(Thing()) == "thing"

    # a real base class wins over a less specific one, like with eager registration:
    # registering the lazy type as a virtual base would put it behind `Mapping`
    Mapping.register(Thing)
    f._clear_cache()
    assert f(Thing()) == "thing"
