## Note on version requirements

The regression test initially passed in CI due to version dependencies. **This bug only manifests with pandas 3.0+ AND numpy < 2.2.3**.

### Why pandas 3.0 matters

anndata has been actively working to support pandas 3.0 (see #2133: "Support pandas 3.0 upcoming changes"). The bug affects users on pandas 3.0 with numpy < 2.2.3.

### Version matrix

| pandas | numpy   | Bug manifests? |
|--------|---------|----------------|
| 2.2.x  | any     | No             |
| 3.0.x  | < 2.2.3 | **Yes**        |
| 3.0.x  | >= 2.2.3| No             |

### Reproducing the bug

To demonstrate the bug exists on `main`, the test must run with:
- pandas >= 3.0
- numpy < 2.2.3

The workflow uses `hatch-test.min` with pandas 3.0 and numpy < 2.2.3 force-installed to reproduce the bug.
