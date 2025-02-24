from __future__ import annotations

import numpy as np
import pytest

import anndata as ad
from anndata._core import extensions


def test_find_stacklevel():
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
        def __init__(self, instance):
            self.instance = instance

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
    assert ns_instance.instance is dummy_obj

    # __get__ should cache the namespace instance on the object.
    # Subsequent access should return the same cached instance.
    assert getattr(dummy_obj, "dummy") is ns_instance


def test_register_namespace(monkeypatch):
    """Test the behavior of the register_anndata_namespace decorator.

    This test verifies that:
    - A new namespace can be registered successfully.
    - The accessor is available on AnnData instances.
    - The accessor is cached on the AnnData instance.
    - An warning is raised if the namespace name is overridden.
    """

    original_dummy = getattr(ad.AnnData, "dummy", None)

    # Register a new namespace called 'dummy'.
    @extensions.register_anndata_namespace("dummy")
    class DummyNamespace:
        def __init__(self, adata: ad.AnnData):
            self.adata = adata

        def greet(self) -> str:
            return "hello"

    # Create an AnnData instance with minimal data.
    rng = np.random.default_rng(42)
    adata = ad.AnnData(X=rng.poisson(1, size=(10, 10)))

    # The accessor should now be available.
    ns_instance = adata.dummy
    assert ns_instance.adata is adata
    assert ns_instance.greet() == "hello"

    # Verify caching behavior on the AnnData instance.
    assert adata.dummy is ns_instance

    # Now, override the same namespace and check that a warning is emitted.
    with pytest.warns(
        UserWarning, match="Overriding existing custom namespace 'dummy'"
    ):

        @extensions.register_anndata_namespace("dummy")
        class DummyNamespaceOverride:
            def __init__(self, adata: ad.AnnData):
                self.adata = adata

            def greet(self) -> str:
                # Return a different string to confirm the override.
                return "world"

    # A new AnnData instance should now use the overridden accessor.
    adata2 = ad.AnnData(X=rng.poisson(1, size=(10, 10)))
    assert adata2.dummy.greet() == "world"

    # Clean up by restoring any previously existing attribute.
    if original_dummy is not None:
        setattr(ad.AnnData, "dummy", original_dummy)
    else:
        if hasattr(ad.AnnData, "dummy"):
            delattr(ad.AnnData, "dummy")


def test_register_reserved_namespace(monkeypatch):
    """
    Check that attempting to register a namespace with a reserved name
    raises an AttributeError.
    """
    reserved_name = "reserved_namespace"

    # Create a new reserved set that includes our test name.
    new_reserved = extensions._reserved_namespaces.union({reserved_name})
    monkeypatch.setattr(extensions, "_reserved_namespaces", new_reserved)

    with pytest.raises(
        AttributeError, match=f"cannot override reserved namespace {reserved_name!r}"
    ):

        @extensions.register_anndata_namespace(reserved_name)
        class DummyNamespace:
            def __init__(self, adata: ad.AnnData):
                self.adata = adata
