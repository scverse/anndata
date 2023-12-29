from __future__ import annotations

import os
import warnings
from typing import TYPE_CHECKING, Any, NamedTuple, TypeVar

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

T = TypeVar("T")


# Heavily inspired by pandas' options mechanism, but stripped down.
# https://github.com/pandas-dev/pandas/blob/86488f700ead42d75dae169032a010c6058e0602/pandas/_config/config.py


class DeprecatedOption(NamedTuple):
    key: str
    msg: str | None
    removal_ver: str | None


# TODO: inherit from Generic[T] as well after python 3.9 is no longer supported
class RegisteredOption(NamedTuple):
    key: str
    defval: object
    doc: str
    validator: Callable[[T], None] | None


# These objects and _describe_option needs to come before docstring_parameterize so they can be used by it and the imports it decorates.
_registered_options: dict[str, RegisteredOption] = {}
_deprecated_options: dict[str, DeprecatedOption] = {}
config: dict[str, object] = {}


def _describe_option(option: str | None = None, print_description=True) -> str:
    if option is not None:
        doc = _registered_options[option].doc
        if option in _deprecated_options:
            doc += "\n"
            opt = _deprecated_options[option]
            doc += opt.msg + "\n"
            doc += f"{option} will be removed in {opt.removal_ver}"
        if print_description:
            print(doc)
        return doc
    else:
        return "\n".join(
            [_describe_option(k, print_description) for k in _registered_options]
        )


def _register_option(
    key: str, defval: object, doc: str, validator: Callable[[T], None]
):
    """Register an option so it can be set/described etc. by end-users

    Parameters
    ----------
    key : str
        Option to be set.
    defval : object
        Default value with which to set the option.
    doc : str
        Docstring for the option
    validator : Callable[[object], Any]
        A function which asserts that the option's value is valid.
    """
    validator(defval)
    _registered_options[key] = RegisteredOption(key, defval, doc, validator)
    config[key] = defval


def docstring_parameterize():
    """This decorator injects the current options into the docstrings of wrapped functions."""

    def dec(obj):
        options_description = _describe_option(print_description=False)
        obj.__doc__ = obj.__doc__.format(
            options_description=options_description,
            available_options=str(", ".join(list(_registered_options.keys()))),
        )
        return obj

    return dec


@docstring_parameterize()
def set_option(key: str, val: object):
    """
    Set an option to a value.  To see the allowed options to be set and their description,
    use describe_option.

    Available options:

    {available_options}

    Parameters
    ----------
    key
        Option to be set.
    val
        Value with which to set the option.

    Raises
    ------
    KeyError
        If the option has not been registered, this function will raise an error.

    Options descriptions:

    {options_description}
    """
    if key not in _registered_options:
        raise KeyError(
            f"{key} is not an available option for anndata.\
            Please open an issue if you believe this is a mistake."
        )
    option = _registered_options[key]
    option.validator(val)
    config[key] = val


@docstring_parameterize()
def get_option(option: str) -> object:
    """
    Gets the option's value.

    Available options:

    {available_options}

    Parameters
    ----------
    option
        Option to be got.

    Returns
    -------
    object
        Value of the option.

    Options descriptions:

    {options_description}
    """
    return config[option]


@docstring_parameterize()
def reset_option(option: str):
    """
    Resets an option to its default value.

    Available options:

    {available_options}

    Parameters
    ----------
    option
        The option to be reset.

    Options descriptions:

    {options_description}
    """
    config[option] = _registered_options[option].defval


@docstring_parameterize()
def describe_option(option: str | None = None, print_description=True) -> str:
    """
    Describe and print (optional) the option(s).

    Available options:

    {available_options}

    Parameters
    ----------
    option
        Option to be described, by default None
    print_description
        Whether or not to also print the description, by default True

    Returns
    -------
    str
        The description

    Options descriptions:

    {options_description}
    """
    return _describe_option(option, print_description)


def check_and_get_environ_var(
    key: str,
    default_value: str,
    allowed_values: Sequence[str] | None = None,
    cast: Callable[[Any], T] = lambda x: x,
) -> T:
    """Get the environment variable and return it is a (potentially) non-string, usable value.

    Parameters
    ----------
    key : str
        The environment variable name.
    default_value : str
        The default value for `os.environ.get`.
    allowed_values : Sequence[str] | None, optional
        Allowable string values., by default None
    cast : _type_, optional
        Casting from the string to a (potentially different) python object, by default lambdax:x

    Returns
    -------
    object
        The casted value.
    """
    environ_val = os.environ.get(key, default_value)
    if allowed_values is not None and environ_val not in allowed_values:
        warnings.warn(
            f'Value "{environ_val}" is not in allowed {allowed_values} for environment variable {key}.\
                      Default {default_value} will be used.'
        )
        return cast(default_value)
    return cast(environ_val)
