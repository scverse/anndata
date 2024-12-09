# Contributing

AnnData follows the development practices outlined in the [Scanpy contribution guide](https://scanpy.readthedocs.io/en/latest/dev/release.html).

```{eval-rst}
.. include:: _key_contributors.rst
```

## CI

### GPU CI

To test GPU specific code we have a paid self-hosted runner to run the gpu specific tests on.
This CI runs by default on the main branch, but for PRs requires the `run-gpu-ci` label to prevent unnecessary runs.
