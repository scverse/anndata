from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

import anndata as ad
from anndata._core import extensions

if TYPE_CHECKING:
    from collections.abc import Generator


@pytest.fixture(autouse=True)
def _cleanup_dummy() -> Generator[None, None, None]:
    """Automatically cleanup dummy namespace after each test."""
    original = getattr(ad.AnnData, "dummy", None)
    yield
    if original is not None:
        ad.AnnData.dummy = original
    elif hasattr(ad.AnnData, "dummy"):
        delattr(ad.AnnData, "dummy")


@pytest.fixture
def dummy_namespace() -> type:
    """Create a basic dummy namespace class."""
    ad.AnnData._accessors = set()

    @ad.register_anndata_namespace("dummy")
    class DummyNamespace:
        def __init__(self, adata: ad.AnnData) -> None:
            self._adata = adata

        def greet(self) -> str:
            return "hello"

    return DummyNamespace


@pytest.fixture
def adata() -> ad.AnnData:
    """Create a basic AnnData object for testing."""
    rng = np.random.default_rng(42)
    return ad.AnnData(X=rng.poisson(1, size=(10, 10)))


def test_find_stacklevel() -> None:
    """Test that find_stacklevel returns a positive integer.

    This function helps determine the correct stacklevel for warnings, so
    we just need to verify it returns a sensible value.
    """
    level = extensions.find_stacklevel()
    assert isinstance(level, int)
    # It should be at least 1, otherwise something is wrong.
    assert level > 0


def test_accessor_namespace() -> None:
    """Test the behavior of the AccessorNameSpace descriptor.

    This test verifies that:
    - When accessed at the class level (i.e., without an instance), the descriptor
      returns the namespace type.
    - When accessed via an instance, the descriptor instantiates the namespace,
      passing the instance to its constructor.
    - The instantiated namespace is then cached on the instance such that subsequent
      accesses of the same attribute return the cached namespace instance.
    """

    # Define a dummy namespace class to be used via the descriptor.
    class DummyNamespace:
        def __init__(self, adata: ad.AnnData) -> None:
            self._adata = adata

        def foo(self) -> str:
            return "foo"

    class Dummy:
        pass

    descriptor = extensions.AccessorNameSpace("dummy", DummyNamespace)

    # When accessed on the class, it should return the namespace type.
    ns_class = descriptor.__get__(None, Dummy)
    assert ns_class is DummyNamespace

    # When accessed via an instance, it should instantiate DummyNamespace.
    dummy_obj = Dummy()
    ns_instance = descriptor.__get__(dummy_obj, Dummy)
    assert isinstance(ns_instance, DummyNamespace)
    assert ns_instance._adata is dummy_obj

    # __get__ should cache the namespace instance on the object.
    # Subsequent access should return the same cached instance.
    assert dummy_obj.dummy is ns_instance


def test_descriptor_instance_caching(dummy_namespace: type, adata: ad.AnnData) -> None:
    """Test that namespace instances are cached on individual AnnData objects."""
    # First access creates the instance
    ns_instance = adata.dummy
    # Subsequent accesses should return the same instance
    assert adata.dummy is ns_instance


def test_register_namespace_basic(dummy_namespace: type, adata: ad.AnnData) -> None:
    """Test basic namespace registration and access."""
    assert adata.dummy.greet() == "hello"


def test_register_namespace_override(dummy_namespace: type) -> None:
    """Test namespace registration and override behavior."""
    assert "dummy" in ad.AnnData._accessors

    # Override should warn and update the namespace
    with pytest.warns(
        UserWarning, match="Overriding existing custom namespace 'dummy'"
    ):

        @ad.register_anndata_namespace("dummy")
        class DummyNamespaceOverride:
            def __init__(self, adata: ad.AnnData) -> None:
                self._adata = adata

            def greet(self) -> str:
                return "world"

    # Verify the override worked
    adata = ad.AnnData(X=np.random.poisson(1, size=(10, 10)))
    assert adata.dummy.greet() == "world"


@pytest.mark.parametrize(
    "attr",
    ["X", "obs", "var", "uns", "obsm", "varm", "layers", "copy", "write"],
)
def test_register_existing_attributes(attr: str) -> None:
    """
    Test that registering an accessor with a name that is a reserved attribute of AnnData raises an attribute error.

    We only test a representative sample of important attributes rather than all of them.
    """
    # Test a representative sample of key AnnData attributes
    with pytest.raises(
        AttributeError,
        match=f"cannot override reserved attribute {attr!r}",
    ):

        @ad.register_anndata_namespace(attr)
        class DummyNamespace:
            def __init__(self, adata: ad.AnnData) -> None:
                self._adata = adata


def test_valid_signature() -> None:
    """Test that a namespace with valid signature is accepted."""

    @ad.register_anndata_namespace("valid")
    class ValidNamespace:
        def __init__(self, adata: ad.AnnData) -> None:
            self.adata = adata


def test_missing_param() -> None:
    """Test that a namespace missing the second parameter is rejected."""
    with pytest.raises(
        TypeError,
        match="Namespace initializer must accept an AnnData instance as the second parameter.",
    ):

        @ad.register_anndata_namespace("missing_param")
        class MissingParamNamespace:
            def __init__(self) -> None:
                pass


def test_wrong_name() -> None:
    """Test that a namespace with wrong parameter name is rejected."""
    with pytest.raises(
        TypeError,
        match="Namespace initializer's second parameter must be named 'adata', got 'notadata'.",
    ):

        @ad.register_anndata_namespace("wrong_name")
        class WrongNameNamespace:
            def __init__(self, notadata: ad.AnnData) -> None:
                self.notadata = notadata


def test_wrong_annotation() -> None:
    """Test that a namespace with wrong parameter annotation is rejected."""
    with pytest.raises(
        TypeError,
        match="Namespace initializer's second parameter must be annotated as the 'AnnData' class, got 'int'.",
    ):

        @ad.register_anndata_namespace("wrong_annotation")
        class WrongAnnotationNamespace:
            def __init__(self, adata: int) -> None:
                self.adata = adata


def test_missing_annotation() -> None:
    """Test that a namespace with missing parameter annotation is rejected."""
    with pytest.raises(AttributeError):

        @ad.register_anndata_namespace("missing_annotation")
        class MissingAnnotationNamespace:
            def __init__(self, adata) -> None:
                self.adata = adata


def test_both_wrong() -> None:
    """Test that a namespace with both wrong name and annotation is rejected."""
    with pytest.raises(
        TypeError,
        match=(
            r"Namespace initializer's second parameter must be named 'adata', got 'info'\. "
            r"And must be annotated as 'AnnData', got 'str'\."
        ),
    ):

        @ad.register_anndata_namespace("both_wrong")
        class BothWrongNamespace:
            def __init__(self, info: str) -> None:
                self.info = info
