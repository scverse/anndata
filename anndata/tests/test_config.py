from __future__ import annotations

import os

import pytest

from anndata._config import (
    SettingsManager,
    check_and_get_environ_var,
)

test_option_doc = "doc string!"
test_option_env_var = "ANNDATA_TEST_VAR"
test_option = "test_var"
default_val = False
test_doc = f"""\
{test_option}: bool
    My doc string!
"""

test_option_doc_2 = "doc string 2!"
test_option_env_var_2 = "ANNDATA_TEST_VAR 2"
test_option_2 = "test_var_2"
default_val_2 = False
test_doc_2 = f"""\
{test_option_2}: bool
    My doc string 2!
"""


def validate_bool(val, option):
    assert val in [True, False], f"{val} not valid boolean for option {option}."


settings = SettingsManager()
settings.register(
    test_option, default_val, test_doc, lambda v: validate_bool(v, test_option)
)

settings.register(
    test_option_2, default_val_2, test_doc_2, lambda v: validate_bool(v, test_option_2)
)


def test_check_and_get_environ_var():
    assert not check_and_get_environ_var(
        test_option_env_var, str(default_val), {"True", "False"}, lambda x: x == "True"
    )
    os.environ[test_option_env_var] = "True"
    assert check_and_get_environ_var(
        test_option_env_var, str(default_val), {"True", "False"}, lambda x: x == "True"
    )
    os.environ[test_option_env_var] = "Not a bool!"
    with pytest.warns(
        match=f'Value "{os.environ[test_option_env_var]}" is not in allowed'
    ):
        check_and_get_environ_var(
            test_option_env_var,
            str(default_val),
            {"True", "False"},
            lambda x: x == "True",
        )


def test_register_option_default():
    assert getattr(settings, test_option) == default_val
    assert settings.describe(test_option) == test_doc


def test_set_option():
    setattr(settings, test_option, not default_val)
    assert getattr(settings, test_option) == (not default_val)
    settings.reset(test_option)
    assert getattr(settings, test_option) == default_val


def test_reset_multiple():
    setattr(settings, test_option, not default_val)
    setattr(settings, test_option_2, not default_val_2)
    settings.reset([test_option, test_option_2])
    assert getattr(settings, test_option) == default_val
    assert getattr(settings, test_option_2) == default_val_2


def test_get_unregistered_option():
    with pytest.raises(AttributeError):
        setattr(settings, test_option + "_different", default_val)


def test_override():
    with settings.override(**{test_option: not default_val}):
        assert getattr(settings, test_option) == (not default_val)
    assert getattr(settings, test_option) == default_val


def test_override_multiple():
    with settings.override(
        **{test_option: not default_val, test_option_2: not default_val_2}
    ):
        assert getattr(settings, test_option) == (not default_val)
        assert getattr(settings, test_option_2) == (not default_val_2)
    assert getattr(settings, test_option) == default_val
    assert getattr(settings, test_option_2) == default_val_2


def test_deprecation():
    warning = "This is a deprecation warning!"
    version = "0.1.0"
    settings.deprecate(test_option, version, warning)
    described_option = settings.describe(test_option, print_description=False).split(
        "\n"
    )
    # first line is message, second two from deprecation
    assert len(described_option) == 5
    assert described_option[-2] == warning
    assert described_option[-1] == f"{test_option} will be removed in {version}"
    with pytest.raises(DeprecationWarning):
        getattr(settings, test_option)


def test_deprecation_no_message():
    version = "0.1.0"
    settings.deprecate(test_option, version)
    described_option = settings.describe(test_option, print_description=False).split(
        "\n"
    )
    # first line is message, second from deprecation version
    assert len(described_option) == 4
    assert described_option[-1] == f"{test_option} will be removed in {version}"
