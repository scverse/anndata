from __future__ import annotations

import inspect
import os
import sys
import textwrap
import warnings
from collections.abc import Iterable
from contextlib import contextmanager
from dataclasses import dataclass, field, fields
from enum import Enum
from functools import partial
from inspect import Parameter, signature
from types import GenericAlias
from typing import TYPE_CHECKING, Generic, NamedTuple, TypeVar, cast

from anndata.compat import CAN_USE_SPARSE_ARRAY
from anndata.compat.exceptiongroups import add_note

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from typing import Any, TypeGuard

T = TypeVar("T")


class DeprecatedOption(NamedTuple):
    option: str
    message: str | None
    removal_version: str | None


def _is_plain_type(obj: object) -> TypeGuard[type]:
    return isinstance(obj, type) and not isinstance(obj, GenericAlias)


def describe(self: RegisteredOption, *, as_rst: bool = False) -> str:
    type_str = self.type.__name__ if _is_plain_type(self.type) else str(self.type)
    if as_rst:
        default_str = repr(self.default_value).replace("\\", "\\\\")
        doc = f"""\
        .. attribute:: settings.{self.option}
           :type: {type_str}
           :value: {default_str}

           {self.description}
        """
    else:
        doc = f"""\
        {self.option}: `{type_str}`
            {self.description} (default: `{self.default_value!r}`).
        """
    return textwrap.dedent(doc)


if sys.version_info >= (3, 11):

    class RegisteredOption(NamedTuple, Generic[T]):
        option: str
        default_value: T
        description: str
        validate: Callable[[T], None]
        type: object

        describe = describe

else:

    class RegisteredOption(NamedTuple):
        option: str
        default_value: T
        description: str
        validate: Callable[[T], None]
        type: object

        describe = describe


def check_and_get_environ_var(
    key: str,
    default_value: str,
    allowed_values: Sequence[str] | None = None,
    cast: Callable[[Any], T] | type[Enum] = lambda x: x,
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
    The casted value.
    """
    environ_value_or_default_value = os.environ.get(key, default_value)
    if (
        allowed_values is not None
        and environ_value_or_default_value not in allowed_values
    ):
        msg = (
            f"Value {environ_value_or_default_value!r} is not in allowed {allowed_values} for environment variable {key}. "
            f"Default {default_value} will be used."
        )
        warnings.warn(msg)
        environ_value_or_default_value = default_value
    return (
        cast(environ_value_or_default_value)
        if not isinstance(cast, type(Enum))
        else cast[environ_value_or_default_value]
    )


def check_and_get_bool(option, default_value):
    return check_and_get_environ_var(
        f"ANNDATA_{option.upper()}",
        str(int(default_value)),
        ["0", "1"],
        lambda x: bool(int(x)),
    )


_docstring = """
This manager allows users to customize settings for the anndata package.
Settings here will generally be for advanced use-cases and should be used with caution.

The following options are available:

{options_description}

For setting an option please use :func:`~anndata.settings.override` (local) or set the above attributes directly (global) i.e., `anndata.settings.my_setting = foo`.
For assignment by environment variable, use the variable name in all caps with `ANNDATA_` as the prefix before import of :mod:`anndata`.
For boolean environment variable setting, use 1 for `True` and 0 for `False`.
"""


@dataclass
class SettingsManager:
    _registered_options: dict[str, RegisteredOption] = field(default_factory=dict)
    _deprecated_options: dict[str, DeprecatedOption] = field(default_factory=dict)
    _config: dict[str, object] = field(default_factory=dict)
    __doc_tmpl__: str = _docstring

    def describe(
        self,
        option: str | Iterable[str] | None = None,
        *,
        should_print_description: bool = True,
        as_rst: bool = False,
    ) -> str:
        """Print and/or return a (string) description of the option(s).

        Parameters
        ----------
        option
            Option(s) to be described, by default None (i.e., do all option)
        should_print_description
            Whether or not to print the description in addition to returning it.

        Returns
        -------
        The description.
        """
        describe = partial(
            self.describe,
            should_print_description=should_print_description,
            as_rst=as_rst,
        )
        if option is None:
            return describe(self._registered_options.keys())
        if isinstance(option, Iterable) and not isinstance(option, str):
            return "\n".join([describe(k) for k in option])
        registered_option = self._registered_options[option]
        doc = registered_option.describe(as_rst=as_rst).rstrip("\n")
        if option in self._deprecated_options:
            opt = self._deprecated_options[option]
            if opt.message is not None:
                doc += f" *{opt.message}"
            doc += f" {option} will be removed in {opt.removal_version}.*"
        if should_print_description:
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
        default_value: T,
        description: str,
        validate: Callable[[T], None],
        option_type: object | None = None,
        get_from_env: Callable[[str, T], T] = lambda x, y: y,
    ) -> None:
        """Register an option so it can be set/described etc. by end-users

        Parameters
        ----------
        option
            Option to be set.
        default_value
            Default value with which to set the option.
        description
            Description to be used in the docstring.
        validate
            A function which raises a `ValueError` or `TypeError` if the value is invalid.
        option_type
            Optional override for the option type to be displayed.  Otherwise `type(default_value)`.
        get_from_env
            An optional function which takes as arguments the name of the option and a default value and returns the value from the environment variable `ANNDATA_CAPS_OPTION` (or default if not present).
            Default behavior is to return `default_value` without checking the environment.
        """
        try:
            validate(default_value)
        except (ValueError, TypeError) as e:
            add_note(e, f"for option {repr(option)}")
            raise e
        option_type = type(default_value) if option_type is None else option_type
        self._registered_options[option] = RegisteredOption(
            option, default_value, description, validate, option_type
        )
        self._config[option] = get_from_env(option, default_value)
        self._update_override_function_for_new_option(option)

    def _update_override_function_for_new_option(
        self,
        option: str,
    ):
        """This function updates the keyword arguments, docstring, and annotations of the `SettingsManager.override` function as the `SettingsManager.register` method is called.

        Parameters
        ----------
        option
            The option being registered for which the override function needs updating.
        """
        option_type = self._registered_options[option].type
        # Update annotations for type checking.
        self.override.__annotations__[option] = option_type
        # __signature__ needs to be updated for tab autocompletion in IPython.
        # See https://github.com/ipython/ipython/issues/11624 for inspiration.
        self.override.__func__.__signature__ = signature(self.override).replace(
            parameters=[
                Parameter(name="self", kind=Parameter.POSITIONAL_ONLY),
                *[
                    Parameter(
                        name=k,
                        annotation=option_type,
                        kind=Parameter.KEYWORD_ONLY,
                    )
                    for k in self._registered_options
                ],
            ]
        )
        # Update docstring for `SettingsManager.override` as well.
        doc = cast(str, self.override.__doc__)
        insert_index = doc.find("\n        Yields")
        option_docstring = "\t" + "\t".join(
            self.describe(option, should_print_description=False).splitlines(
                keepends=True
            )
        )
        self.override.__func__.__doc__ = (
            f"{doc[:insert_index]}\n{option_docstring}{doc[insert_index:]}"
        )

    def __setattr__(self, option: str, val: object) -> None:
        """
        Set an option to a value.  To see the allowed option to be set and their description,
        use describe_option.

        Parameters
        ----------
        option
            Option to be set.
        val
            Value with which to set the option.

        Raises
        ------
        AttributeError
            If the option has not been registered, this function will raise an error.
        """
        if option in {f.name for f in fields(self)}:
            return super().__setattr__(option, val)
        elif option not in self._registered_options:
            msg = (
                f"{option} is not an available option for anndata. "
                "Please open an issue if you believe this is a mistake."
            )
            raise AttributeError(msg)
        registered_option = self._registered_options[option]
        registered_option.validate(val)
        self._config[option] = val

    def __getattr__(self, option: str) -> object:
        """
        Gets the option's value.

        Parameters
        ----------
        option
            Option to be got.

        Returns
        -------
        Value of the option.
        """
        if option in self._deprecated_options:
            deprecated = self._deprecated_options[option]
            msg = f"{repr(option)} will be removed in {deprecated.removal_version}. {deprecated.message}"
            warnings.warn(msg, DeprecationWarning)
        if option in self._config:
            return self._config[option]
        msg = f"{option} not found."
        raise AttributeError(msg)

    def __dir__(self) -> Iterable[str]:
        return sorted((*dir(super()), *self._config.keys()))

    def reset(self, option: Iterable[str] | str) -> None:
        """
        Resets option(s) to its (their) default value(s).

        Parameters
        ----------
        option
            The option(s) to be reset.
        """
        if isinstance(option, Iterable) and not isinstance(option, str):
            for opt in option:
                self.reset(opt)
        else:
            self._config[option] = self._registered_options[option].default_value

    @contextmanager
    def override(self, **overrides):
        """
        Provides local override via keyword arguments as a context manager.

        Parameters
        ----------

        Yields
        ------
        None
        """
        restore = {a: getattr(self, a) for a in overrides}
        try:
            for attr, value in overrides.items():
                setattr(self, attr, value)
            yield None
        finally:
            for attr, value in restore.items():
                setattr(self, attr, value)

    def __repr__(self) -> str:
        params = "".join(f"\t{k}={v!r},\n" for k, v in self._config.items())
        return f"{type(self).__name__}(\n{params}\n)"

    @property
    def __doc__(self):
        in_sphinx = any("/sphinx/" in frame.filename for frame in inspect.stack())
        options_description = self.describe(
            should_print_description=False, as_rst=in_sphinx
        )
        return self.__doc_tmpl__.format(
            options_description=options_description,
        )


settings = SettingsManager()

##################################################################################
# PLACE REGISTERED SETTINGS HERE SO THEY CAN BE PICKED UP FOR DOCSTRING CREATION #
##################################################################################


def validate_bool(val: Any) -> None:
    if not isinstance(val, bool):
        msg = f"{val} not valid boolean"
        raise TypeError(msg)


settings.register(
    "remove_unused_categories",
    default_value=True,
    description="Whether or not to remove unused categories with :class:`~pandas.Categorical`.",
    validate=validate_bool,
    get_from_env=check_and_get_bool,
)

settings.register(
    "check_uniqueness",
    default_value=True,
    description=(
        "Whether or not to check uniqueness of the `obs` indices on `__init__` of :class:`~anndata.AnnData`."
    ),
    validate=validate_bool,
    get_from_env=check_and_get_bool,
)

settings.register(
    "allow_write_nullable_strings",
    default_value=False,
    description="Whether or not to allow writing of `pd.arrays.StringArray`.",
    validate=validate_bool,
    get_from_env=check_and_get_bool,
)


def validate_sparse_settings(val: Any) -> None:
    validate_bool(val)
    if not CAN_USE_SPARSE_ARRAY and cast(bool, val):
        msg = (
            "scipy.sparse.cs{r,c}array is not available in current scipy version. "
            "Falling back to scipy.sparse.cs{r,c}_matrix for reading."
        )
        raise ValueError(msg)


settings.register(
    "use_sparse_array_on_read",
    default_value=False,
    description="Whether or not to use :class:`scipy.sparse.sparray` as the default class when reading in data",
    validate=validate_sparse_settings,
    get_from_env=check_and_get_bool,
)

##################################################################################
##################################################################################
