# Contributing

AnnData follows the development practices outlined in the [Scanpy contribution guide](https://scanpy.readthedocs.io/en/latest/dev/release.html).

```{eval-rst}
.. include:: _key_contributors.rst
```

## Code checks

We use [prek][] to run the [pre-commit][] hooks that enforce code style, formatting, and static type checking.
Install it and run `prek install` once to set up the git hook, then `prek run --all-files` to check the whole repository.
The same hooks run in the `Pre-commit checks` job of the CI, and Dependabot keeps their versions current.

Static type checking is done by [mypy][] over `src` and `tests`, configured through `[tool.mypy]` in the `pyproject.toml`.
Type hints are optional: mypy only checks functions that carry annotations, so untyped code can be annotated gradually.

[prek]: https://prek.j178.dev/
[pre-commit]: https://pre-commit.com/
[mypy]: https://mypy.readthedocs.io/

## CI

### GPU CI

To test GPU specific code we have a paid self-hosted runner to run the gpu specific tests on.
This CI runs by default on the main branch, but for PRs requires the `run-gpu-ci` label to prevent unnecessary runs.
