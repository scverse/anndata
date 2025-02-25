from __future__ import annotations

import numpy as np
import pytest

import anndata as ad
from anndata._core import extensions


def test_find_stacklevel():
    """Test that find_stacklevel returns a positive integer.

    This function helps determine the correct stacklevel for warnings, so
    we just need to verify it returns a sensible value.
    """
    level = extensions.find_stacklevel()
    assert isinstance(level, int)
    # It should be at least 1, otherwise something is wrong.
    assert level > 0


def test_accessor_namespace():
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
        def __init__(self, adata: ad.AnnData):
            self._adata = adata

        def foo(self):
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
    assert getattr(dummy_obj, "dummy") is ns_instance


def test_register_namespace_basic():
    """Test the basic behavior of the register_anndata_namespace decorator.

    This test verifies that:
    - A new namespace can be registered successfully.
    - The accessor is available on AnnData instances.
    """
    original_dummy = getattr(ad.AnnData, "dummy", None)

    # Register a new namespace called 'dummy'.
    @ad.register_anndata_namespace("dummy")
    class DummyNamespace:
        def __init__(self, adata: ad.AnnData):
            self._adata = adata

        def greet(self) -> str:
            return "hello"

    # Create an AnnData instance with minimal data.
    rng = np.random.default_rng(42)
    adata = ad.AnnData(X=rng.poisson(1, size=(10, 10)))

    # The accessor should now be available.
    ns_instance = adata.dummy
    assert ns_instance._adata is adata
    assert ns_instance.greet() == "hello"

    # Clean up
    if original_dummy is not None:
        setattr(ad.AnnData, "dummy", original_dummy)
    else:
        if hasattr(ad.AnnData, "dummy"):
            delattr(ad.AnnData, "dummy")


def test_register_namespace_caching():
    """Test that the namespace accessor is cached on the AnnData instance."""
    original_dummy = getattr(ad.AnnData, "dummy", None)

    # Register a new namespace
    @ad.register_anndata_namespace("dummy")
    class DummyNamespace:
        def __init__(self, adata: ad.AnnData):
            self._adata = adata

        def greet(self) -> str:
            return "hello"

    # Create an AnnData instance
    rng = np.random.default_rng(42)
    adata = ad.AnnData(X=rng.poisson(1, size=(10, 10)))

    # Access the namespace to trigger caching
    ns_instance = adata.dummy

    # Verify caching behavior on the AnnData instance.
    assert adata.dummy is ns_instance

    # Clean up
    if original_dummy is not None:
        setattr(ad.AnnData, "dummy", original_dummy)
    else:
        if hasattr(ad.AnnData, "dummy"):
            delattr(ad.AnnData, "dummy")


def test_register_namespace_override():
    """Test that a warning is raised when overriding an existing namespace."""
    original_dummy = getattr(ad.AnnData, "dummy", None)

    # Register a namespace first
    @ad.register_anndata_namespace("dummy")
    class DummyNamespace:
        def __init__(self, adata: ad.AnnData):
            self._adata = adata

        def greet(self) -> str:
            return "hello"

    # Now, override the same namespace and check that a warning is emitted.
    with pytest.warns(
        UserWarning, match="Overriding existing custom namespace 'dummy'"
    ):

        @ad.register_anndata_namespace("dummy")
        class DummyNamespaceOverride:
            def __init__(self, adata: ad.AnnData):
                self._adata = adata

            def greet(self) -> str:
                # Return a different string to confirm the override.
                return "world"

    # A new AnnData instance should now use the overridden accessor.
    rng = np.random.default_rng(42)
    adata = ad.AnnData(X=rng.poisson(1, size=(10, 10)))
    assert adata.dummy.greet() == "world"

    # Clean up
    if original_dummy is not None:
        setattr(ad.AnnData, "dummy", original_dummy)
    else:
        if hasattr(ad.AnnData, "dummy"):
            delattr(ad.AnnData, "dummy")


def test_register_existing_attributes():
    """
    Test that registering an accessor with a name that is a reserved attribute of AnnData raises an attribute error.

    We only test a representative sample of important attributes rather than all of them.
    """
    # Test a representative sample of key AnnData attributes
    key_attributes = [
        "X",
        "obs",
        "var",
        "uns",
        "obsm",
        "varm",
        "layers",
        "copy",
        "write",
    ]

    for attr in key_attributes:
        with pytest.raises(
            AttributeError,
            match=f"cannot override reserved attribute {attr!r}",
        ):

            @ad.register_anndata_namespace(attr)
            class DummyNamespace:
                def __init__(self, adata: ad.AnnData) -> None:
                    self._adata = adata


def test_check_namespace_signature_valid():
    """Test that a namespace with valid signature is accepted."""

    # Valid namespace: correct signature.
    # should not raise any error.
    @ad.register_anndata_namespace("valid")
    class ValidNamespace:
        def __init__(self, adata: ad.AnnData) -> None:
            self.adata = adata


def test_check_namespace_signature_missing_param():
    """Test that a namespace missing the second parameter is rejected."""
    with pytest.raises(
        TypeError,
        match="Namespace initializer must accept an AnnData instance as the second parameter.",
    ):

        @ad.register_anndata_namespace("missing_param")
        class MissingParamNamespace:
            def __init__(self) -> None:
                pass


def test_check_namespace_signature_wrong_name():
    """Test that a namespace with wrong parameter name is rejected."""
    with pytest.raises(
        TypeError,
        match="Namespace initializer's second parameter must be named 'adata', got 'notadata'.",
    ):

        @ad.register_anndata_namespace("wrong_name")
        class WrongNameNamespace:
            def __init__(self, notadata: ad.AnnData) -> None:
                self.notadata = notadata


def test_check_namespace_signature_wrong_annotation():
    """Test that a namespace with wrong parameter annotation is rejected."""
    with pytest.raises(
        TypeError,
        match="Namespace initializer's second parameter must be annotated as the 'AnnData' class, got 'int'.",
    ):

        @ad.register_anndata_namespace("wrong_annotation")
        class WrongAnnotationNamespace:
            def __init__(self, adata: int) -> None:
                self.adata = adata


def test_check_namespace_signature_missing_annotation():
    """Test that a namespace with missing parameter annotation is rejected."""
    with pytest.raises(AttributeError):

        @ad.register_anndata_namespace("missing_annotation")
        class MissingAnnotationNamespace:
            def __init__(self, adata) -> None:
                self.adata = adata


def test_check_namespace_signature_both_wrong():
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
