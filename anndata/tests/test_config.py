from __future__ import annotations

import os
from enum import Enum

import pytest

from anndata._settings import (
    SettingsManager,
    check_and_get_bool,
    check_and_get_environ_var,
    validate_bool,
)

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


def test_register_with_env(monkeypatch):
    with monkeypatch.context() as mp:
        option_env = "test_var_env"
        default_val_env = False
        description_env = "My doc string env!"
        option_env_var = "ANNDATA_" + option_env.upper()
        mp.setenv(option_env_var, "1")

        settings.register(
            option_env,
            default_val_env,
            description_env,
            validate_bool,
            get_from_env=check_and_get_bool,
        )

        assert settings.test_var_env


def test_register_with_env_enum(monkeypatch):
    with monkeypatch.context() as mp:
        option_env = "test_var_env"
        default_val_env = False
        description_env = "My doc string env!"
        option_env_var = "ANNDATA_" + option_env.upper()
        mp.setenv(option_env_var, "b")

        class TestEnum(Enum):
            a = False
            b = True

        def check_and_get_bool_enum(option, default_value):
            return check_and_get_environ_var(
                "ANNDATA_" + option.upper(), "a", cast=TestEnum
            ).value

        settings.register(
            option_env,
            default_val_env,
            description_env,
            validate_bool,
            get_from_env=check_and_get_bool_enum,
        )

        assert settings.test_var_env


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


def test_check_and_get_environ_var(monkeypatch):
    with monkeypatch.context() as mp:
        option_env_var = "ANNDATA_OPTION"
        assert hash("foo") == check_and_get_environ_var(
            option_env_var, "foo", ["foo", "bar"], lambda x: hash(x)
        )
        mp.setenv(option_env_var, "bar")
        assert hash("bar") == check_and_get_environ_var(
            option_env_var, "foo", ["foo", "bar"], lambda x: hash(x)
        )
        mp.setenv(option_env_var, "Not foo or bar")
        with pytest.warns(
            match=f'Value "{os.environ[option_env_var]}" is not in allowed'
        ):
            check_and_get_environ_var(
                option_env_var, "foo", ["foo", "bar"], lambda x: hash(x)
            )
        assert hash("Not foo or bar") == check_and_get_environ_var(
            option_env_var, "foo", cast=lambda x: hash(x)
        )


def test_check_and_get_bool(monkeypatch):
    with monkeypatch.context() as mp:
        option_env_var = "ANNDATA_" + option.upper()
        assert not check_and_get_bool(option, default_val)
        mp.setenv(option_env_var, "1")
        assert check_and_get_bool(option, default_val)
        mp.setenv(option_env_var, "Not 0 or 1")
        with pytest.warns(
            match=f'Value "{os.environ[option_env_var]}" is not in allowed'
        ):
            check_and_get_bool(option, default_val)


def test_check_and_get_bool_enum(monkeypatch):
    with monkeypatch.context() as mp:
        option_env_var = "ANNDATA_" + option.upper()
        mp.setenv(option_env_var, "b")

        class TestEnum(Enum):
            a = False
            b = True

        assert check_and_get_environ_var(option_env_var, "a", cast=TestEnum).value
