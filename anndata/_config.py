from __future__ import annotations

import os
import warnings
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any, NamedTuple, TypeVar

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

T = TypeVar("T")


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


class SettingsManager:
    """
    This manager allows user to customize settings for the anndata package.

    Available options:

    {available_options}

    Options descriptions:

    {options_description}

    """

    _registered_options: dict[str, RegisteredOption] = {}
    _deprecated_options: dict[str, DeprecatedOption] = {}
    _config: dict[str, object] = {}

    def describe(self, option: str | None = None, print_description=True) -> str:
        if option is not None:
            doc = self._registered_options[option].doc
            if option in self._deprecated_options:
                doc += "\n"
                opt = self._deprecated_options[option]
                doc += opt.msg + "\n"
                doc += f"{option} will be removed in {opt.removal_ver}"
            if print_description:
                print(doc)
            return doc
        else:
            return "\n".join(
                [self.describe(k, print_description) for k in self._registered_options]
            )

    def register(
        self, key: str, defval: object, doc: str, validator: Callable[[T], None]
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
        self._registered_options[key] = RegisteredOption(key, defval, doc, validator)
        self._config[key] = defval

    def __setitem__(self, key: str, val: object):
        """
        Set an option to a value.  To see the allowed options to be set and their description,
        use describe_option.

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
        """
        if key not in self._registered_options:
            raise KeyError(
                f"{key} is not an available option for anndata.\
                Please open an issue if you believe this is a mistake."
            )
        option = self._registered_options[key]
        option.validator(val)
        self._config[key] = val

    def __getitem__(self, option: str) -> object:
        """
        Gets the option's value.

        Parameters
        ----------
        option
            Option to be got.

        Returns
        -------
        object
            Value of the option.
        """
        return self._config[option]

    def reset(self, option: str):
        """
        Resets an option to its default value.

        Parameters
        ----------
        option
            The option to be reset.
        """
        self._config[option] = self._registered_options[option].defval

    @contextmanager
    def override(self, **overrides):
        """
        Provides local override via keyword arguments.

        Yields
        ------
        None
        """
        restore = {a: self[a] for a in overrides}
        try:
            for attr, value in overrides.items():
                self[attr] = value
            yield None
        finally:
            for attr, value in restore.items():
                self[attr] = value


settings = SettingsManager()

##################################################################################
# PLACE REGISTERED SETTINGS HERE SO THEY CAN BE PICKED UP FOR DOCSTRING CREATION #
##################################################################################


##################################################################################
##################################################################################


options_description = settings.describe(print_description=False)
settings.__doc__ = settings.__doc__.format(
    options_description=options_description,
    available_options=str(", ".join(list(settings._registered_options.keys()))),
)


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
