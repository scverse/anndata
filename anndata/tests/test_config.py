from __future__ import annotations

import os

import pytest

from anndata._config import (
    SettingsManager,
    check_and_get_bool,
    check_and_get_environ_var,
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


def validate_bool(val, option):
    assert val in [True, False], f"{val} not valid boolean for option {option}."


def validate_int_list(val, option):
    assert [
        isinstance(type(e), int) for e in val
    ], f"{val} not valid int list for option {option}."


settings = SettingsManager()
settings.register(option, default_val, description, lambda v: validate_bool(v, option))

settings.register(
    option_2, default_val_2, description_2, lambda v: validate_bool(v, option_2)
)

settings.register(
    option_3,
    default_val_3,
    description_3,
    lambda v: validate_int_list(v, option_3),
    type_3,
)


def test_register_option_default():
    assert getattr(settings, option) == default_val
    assert description in settings.describe(option)


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
    with pytest.raises(DeprecationWarning):
        getattr(settings, option)


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
