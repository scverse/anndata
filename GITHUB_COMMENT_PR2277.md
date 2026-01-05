## Note on version requirements

The regression test initially passed in CI due to version dependencies. **This bug only manifests with pandas 3.0+ AND numpy < 2.2.3**.

### Why these versions matter

- **pandas 3.0**: anndata has been actively working to support pandas 3.0 (see #2133: "Support pandas 3.0 upcoming changes").
- **numpy < 2.2.3**: NumPy 2.2.3 fixed incorrect bytes-to-string coercion ([numpy#28282](https://github.com/numpy/numpy/pull/28282)), which accidentally masks this bug. Before 2.2.3, bytes like `b'cell_A'` were incorrectly converted to `"b'cell_A'"` instead of `"cell_A"`.

### Version matrix

| pandas | numpy   | Bug manifests? |
|--------|---------|----------------|
| 2.2.x  | any     | No             |
| 3.0.x  | < 2.2.3 | **Yes**        |
| 3.0.x  | >= 2.2.3| No (masked by numpy fix) |

### Why the fix in #2272 is still needed

Even though numpy 2.2.3+ masks this bug, the fix in #2272 may still be the correct solution because:
1. Users may still be on numpy < 2.2.3
2. The anndata code should properly decode bytes rather than relying on numpy's coercion behavior
3. Using `.asstr()` explicitly handles HDF5 byte strings correctly

### Reproducing the bug

To demonstrate the bug exists on `main`, the test must run with:
- pandas >= 3.0
- numpy < 2.2.3

The workflow uses `hatch-test.min` with pandas 3.0 and numpy < 2.2.3 force-installed to reproduce the bug.
