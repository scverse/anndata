from __future__ import annotations

import pytest

from anndata._config import SettingsManager

option = "test_var"
default_val = False
description = "My doc string!"

option_2 = "test_var_2"
default_val_2 = False
description_2 = "My doc string 2!"

option_3 = "test_var_3"
default_val_3 = [1, 2]
description_3 = "My doc string 3!"
type_3 = list[int]


def validate_bool(val) -> bool:
    if not isinstance(val, bool):
        raise TypeError(f"{val} not valid boolean")
    return True


def validate_int_list(val) -> bool:
    if not isinstance(val, list) or not [isinstance(type(e), int) for e in val]:
        raise TypeError(f"{repr(val)} is not a valid int list")
    return True


settings = SettingsManager()
settings.register(option, default_val, description, validate_bool)

settings.register(option_2, default_val_2, description_2, validate_bool)

settings.register(
    option_3,
    default_val_3,
    description_3,
    validate_int_list,
    type_3,
)


def test_register_option_default():
    assert getattr(settings, option) == default_val
    assert description in settings.describe(option)


def test_register_bad_option():
    with pytest.raises(TypeError, match="'foo' is not a valid int list"):
        settings.register(
            "test_var_4",
            "foo",  # should be a list of ints
            description_3,
            validate_int_list,
            type_3,
        )


def test_set_option():
    setattr(settings, option, not default_val)
    assert getattr(settings, option) == (not default_val)
    settings.reset(option)
    assert getattr(settings, option) == default_val


def test_dir():
    assert {option, option_2, option_3} <= set(dir(settings))
    assert dir(settings) == sorted(dir(settings))


def test_reset_multiple():
    setattr(settings, option, not default_val)
    setattr(settings, option_2, not default_val_2)
    settings.reset([option, option_2])
    assert getattr(settings, option) == default_val
    assert getattr(settings, option_2) == default_val_2


def test_get_unregistered_option():
    with pytest.raises(AttributeError):
        setattr(settings, option + "_different", default_val)


def test_override():
    with settings.override(**{option: not default_val}):
        assert getattr(settings, option) == (not default_val)
    assert getattr(settings, option) == default_val


def test_override_multiple():
    with settings.override(**{option: not default_val, option_2: not default_val_2}):
        assert getattr(settings, option) == (not default_val)
        assert getattr(settings, option_2) == (not default_val_2)
    assert getattr(settings, option) == default_val
    assert getattr(settings, option_2) == default_val_2


def test_deprecation():
    warning = "This is a deprecation warning!"
    version = "0.1.0"
    settings.deprecate(option, version, warning)
    described_option = settings.describe(option, print_description=False)
    # first line is message, second two from deprecation
    default_deprecation_message = f"{option} will be removed in {version}.*"
    assert described_option.endswith(default_deprecation_message)
    described_option = (
        described_option.rstrip().removesuffix(default_deprecation_message).rstrip()
    )
    assert described_option.endswith(warning)
    with pytest.warns(
        DeprecationWarning,
        match="'test_var' will be removed in 0.1.0. This is a deprecation warning!",
    ):
        assert getattr(settings, option) == default_val


def test_deprecation_no_message():
    version = "0.1.0"
    settings.deprecate(option, version)
    described_option = settings.describe(option, print_description=False)
    # first line is message, second from deprecation version
    assert described_option.endswith(f"{option} will be removed in {version}.*")


def test_option_typing():
    assert settings._registered_options[option_3].type == type_3
    assert str(type_3) in settings.describe(option_3, print_description=False)
