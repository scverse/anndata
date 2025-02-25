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

    # Verify caching behavior on the AnnData instance.
    assert adata.dummy is ns_instance

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
    adata2 = ad.AnnData(X=rng.poisson(1, size=(10, 10)))
    assert adata2.dummy.greet() == "world"

    # Clean up by restoring any previously existing attribute.
    if original_dummy is not None:
        setattr(ad.AnnData, "dummy", original_dummy)
    else:
        if hasattr(ad.AnnData, "dummy"):
            delattr(ad.AnnData, "dummy")


def test_register_existing_attributes(monkeypatch):
    """
    Test that registering an accessor with a name that is an attribute of AnnData raises an attribute error.

    i.e. we do not want users to override say, `AnnData.X` or `AnnData.obs_names`, etc...
    """
    reserved_names = dir(ad.AnnData)

    for reserved_name in reserved_names:
        with pytest.raises(
            AttributeError,
            match=f"cannot override reserved attribute {reserved_name!r}",
        ):

            @ad.register_anndata_namespace(reserved_name)
            class DummyNamespace:
                def __init__(self, adata: ad.AnnData) -> None:
                    self._adata = adata


def test_check_namespace_signature_comprehensive():
    """Comprehensive test for _check_namespace_signature covering all edge cases.

    We test the following cases:
    1. Valid namespace: correct signature.
    2. Missing the second parameter (i.e. the class only has self as a parameter).
    3. Wrong parameter name: second parameter not named 'adata'.
    4. Wrong annotation: second parameter annotated as wrong type.
    5. Both wrong: wrong name and wrong annotation.
    6. Missing annotation: no type annotation provided on the second parameter.

    """

    # 1. Valid namespace: correct signature.
    class ValidNamespace:
        def __init__(self, adata: ad.AnnData) -> None:
            self.adata = adata

    # Should not raise any error.
    extensions._check_namespace_signature(ValidNamespace)

    # 2. Missing the second parameter.
    class MissingParamNamespace:
        def __init__(self) -> None:
            pass

    with pytest.raises(
        TypeError,
        match="Namespace initializer must accept an AnnData instance as the second parameter.",
    ):
        extensions._check_namespace_signature(MissingParamNamespace)

    # 3. Wrong parameter name: second parameter not named 'adata'.
    class WrongNameNamespace:
        def __init__(self, notadata: ad.AnnData) -> None:
            self.notadata = notadata

    with pytest.raises(
        TypeError,
        match="Namespace initializer's second parameter must be named 'adata', got 'notadata'.",
    ):
        extensions._check_namespace_signature(WrongNameNamespace)

    # 4. Wrong annotation: second parameter annotated as wrong type.
    class WrongAnnotationNamespace:
        def __init__(self, adata: int) -> None:
            self.adata = adata

    with pytest.raises(
        TypeError,
        match="Namespace initializer's second parameter must be annotated as the 'AnnData' class, got 'int'.",
    ):
        extensions._check_namespace_signature(WrongAnnotationNamespace)

    # 5. Both wrong: wrong name and wrong annotation.
    class BothWrongNamespace:
        def __init__(self, info: str) -> None:
            self.info = info

    with pytest.raises(
        TypeError,
        match=(
            r"Namespace initializer's second parameter must be named 'adata', got 'info'\. "
            r"And must be annotated as 'AnnData', got 'str'\."
        ),
    ):
        extensions._check_namespace_signature(BothWrongNamespace)

    # 6. Missing annotation: no type annotation provided on the second parameter.
    class MissingAnnotationNamespace:
        def __init__(self, adata) -> None:
            self.adata = adata

    with pytest.raises(AttributeError):
        extensions._check_namespace_signature(MissingAnnotationNamespace)
