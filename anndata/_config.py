from __future__ import annotations

import os
import warnings
from typing import TYPE_CHECKING, Any, Callable, NamedTuple

if TYPE_CHECKING:
    from collections.abc import Sequence

# Heavily inspired by pandas' options mechanism, but stripped down.
# https://github.com/pandas-dev/pandas/blob/86488f700ead42d75dae169032a010c6058e0602/pandas/_config/config.py

config = {}


class DeprecatedOption(NamedTuple):
    key: str
    msg: str | None
    removal_ver: str | None


class RegisteredOption(NamedTuple):
    key: str
    defval: object
    doc: str
    validator: Callable[[object], Any] | None


registered_options: dict[str, RegisteredOption] = {}
deprecated_options: dict[str, DeprecatedOption] = {}
config: dict[str, object] = {}


def register_option(
    key: str, defval: object, doc: str, validator: Callable[[object], Any]
):
    validator(defval)
    registered_options[key] = RegisteredOption(key, defval, doc, validator)
    config[key] = defval


def set_option(key: str, val: object):
    option = registered_options[key]
    option.validator(val)
    config[key] = val


def get_option(option: str) -> object:
    return config[option]


def reset_option(option: str):
    config[option] = registered_options[option].defval


def describe_option(option: str | None = None):
    if option is not None:
        print(registered_options[option].doc)
    else:
        for k in registered_options:
            print(registered_options[k].doc + "\n")
            if k in deprecated_options:
                opt = deprecated_options[k]
                print(opt.msg + "\n")
                print(f"{k} will be removed in {opt.removal_ver}" + "\n")


def check_and_get_environ_var(
    key: str, defaultval: str, options: Sequence[str], cast: Callable[[object], object]
) -> object:
    environ_val = os.environ.get(key, defaultval)
    if environ_val not in options:
        warnings.warn(
            f'Value "{environ_val}" is not in allowed {options} for environment variable {key}.\
                      Default {defaultval} will be used.'
        )
        return cast(defaultval)
    return cast(environ_val)
