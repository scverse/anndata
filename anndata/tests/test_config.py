from __future__ import annotations

import os

import pytest

from anndata._config import (
    _register_option,
    check_and_get_environ_var,
    describe_option,
    get_option,
    reset_option,
    set_option,
)

test_option_doc = "doc string!"
test_option_env_var = "ANNDATA_TEST_VAR"
test_option = "test_var"
default_val = False
test_doc = """My doc string!"""


def validate_bool(val, option):
    assert val in [True, False], f"{val} not valid boolean for option {option}."


_register_option(
    test_option, default_val, test_doc, lambda v: validate_bool(v, test_option)
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


def test__register_option_default():
    assert get_option(test_option) == default_val
    assert describe_option(test_option) == test_doc


def test_set_option():
    set_option(test_option, not default_val)
    assert get_option(test_option) == (not default_val)
    reset_option(test_option)
    assert get_option(test_option) == default_val


def test_get_unregistered_option():
    with pytest.raises(KeyError):
        set_option(test_option + "_different", default_val)
