# Contributing

AnnData follows the development practices outlined in the [Scanpy contribution guide](https://scanpy.readthedocs.io/en/latest/dev/index.html).

```{eval-rst}
.. include:: _key_contributors.rst
```

## Release Notes

AnnData differs from `scanpy` (for now) in how its releases are done.
It uses [towncrier][] to build its changelog.
We have set up some automation around this process.
To run `towncrier`, create a `PR` into the base branch of the release with the compiled changelog, and backport to `main` if needed (i.e., the base branch is something like `0.10.x`), run

```shell
hatch run towncrier:build --version="X.Y.Z"
```

You may add the option `--dry-run` at the end to do the local steps without pushing to Github, although the push will be mocked via [`gh pr --dry-run`](https://cli.github.com/manual/gh_pr_create).

[towncrier]: https://towncrier.readthedocs.io/en/stable/

## CI

### GPU CI

To test GPU specific code we have a paid self-hosted runner to run the gpu specific tests on.
This CI runs by default on the main branch, but for PRs requires the `run-gpu-ci` label to prevent unnecessary runs.
