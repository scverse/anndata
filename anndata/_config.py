from __future__ import annotations

import os
import warnings
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any, NamedTuple, TypeVar

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

T = TypeVar("T")


class DeprecatedOption(NamedTuple):
    option: str
    message: str | None
    removal_version: str | None


# TODO: inherit from Generic[T] as well after python 3.9 is no longer supported
class RegisteredOption(NamedTuple):
    option: str
    default_value: object
    doc: str
    validator: Callable[[T], None] | None


class SettingsManager:
    """
    This manager allows users to customize settings for the anndata package.

    Available options:

    {available_options}

    Options descriptions:

    {options_description}

    """

    _registered_options: dict[str, RegisteredOption] = {}
    _deprecated_options: dict[str, DeprecatedOption] = {}
    _config: dict[str, object] = {}

    def describe(
        self, option: str | None = None, print_description: bool = True
    ) -> str:
        """Print and/or return a (string) description of the option(s).

        Parameters
        ----------
        option
            Option to be described, by default None (i.e., do all options)
        print_description
            Whether or not to print the description in addition to returning it., by default True

        Returns
        -------
        str
            The description.
        """
        if option is None:
            return "\n".join(
                [self.describe(k, print_description) for k in self._registered_options]
            )

        doc = self._registered_options[option].doc
        if option in self._deprecated_options:
            doc += "\n"
            opt = self._deprecated_options[option]
            if opt.message is not None:
                doc += opt.message + "\n"
            doc += f"{option} will be removed in {opt.removal_version}"
        if print_description:
            print(doc)
        return doc

    def deprecate(
        self, option: str, removal_version: str, message: str | None = None
    ) -> None:
        """Deprecate options with a message at a version.

        Parameters
        ----------
        option
            Which option should be deprecated.
        removal_version
            The version targeted for removal.
        message
            A custom message.
        """
        self._deprecated_options[option] = DeprecatedOption(
            option, message, removal_version
        )

    def register(
        self,
        option: str,
        default_value: object,
        doc: str,
        validator: Callable[[T], None],
    ) -> None:
        """Register an option so it can be set/described etc. by end-users

        Parameters
        ----------
        option
            Option to be set.
        default_value
            Default value with which to set the option.
        doc
            Docstring for the option
        validator
            A function which asserts that the option's value is valid.
        """
        validator(default_value)
        self._registered_options[option] = RegisteredOption(
            option, default_value, doc, validator
        )
        self._config[option] = default_value

    def __setitem__(self, option: str, val: object) -> None:
        """
        Set an option to a value.  To see the allowed options to be set and their description,
        use describe_option.

        Parameters
        ----------
        option
            Option to be set.
        val
            Value with which to set the option.

        Raises
        ------
        KeyError
            If the option has not been registered, this function will raise an error.
        """
        if option not in self._registered_options:
            raise KeyError(
                f"{option} is not an available option for anndata.\
                Please open an issue if you believe this is a mistake."
            )
        registered_option = self._registered_options[option]
        registered_option.validator(val)
        self._config[option] = val

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

    def reset(self, option: str) -> None:
        """
        Resets an option to its default value.

        Parameters
        ----------
        option
            The option to be reset.
        """
        self._config[option] = self._registered_options[option].default_value

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
    key
        The environment variable name.
    default_value
        The default value for `os.environ.get`.
    allowed_values
        Allowable string values., by default None
    cast
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
