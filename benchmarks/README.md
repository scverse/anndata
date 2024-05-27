# AnnData Benchmarks

This repo contains some work in progress benchmarks for [AnnData](https://github.com/theislab/anndata) using [asv](https://asv.readthedocs.io).

## Setup

I definitely recommend reading through the asv docs. Currently, this assumes the benchmark suite can reach the `anndata` repo via the path `../anndata`. Otherwise, all you'll need to do is create a [machine file](https://asv.readthedocs.io/en/stable/commands.html#asv-machine) for your system and make sure `anndata`s dependencies are installable via `conda`.

### Data

Data will need to be retrieved for these benchmarks. This can be downloaded using the script fetch_datasets.py.

Note that the `h5ad` format has changed since it's inception. While the `anndata` package maintains backwards compatibility, older versions of `anndata` will not be able to read files written by more recent versions. To get around this for the benchmarks, datasets have to be able to be read by all versions which can require a setup function that creates the anndata object.

## Usage

### Runnings the benchmarks:

To run benchmarks for a particular commit: `asv run {commit} --steps 1 -b`

To run benchmarks for a range of commits: `asv run {commit1}..{commit2}`

You can filter out the benchmarks which are run with the `-b {pattern}` flag.

### Accessing the benchmarks

You can see what benchmarks you've already run using `asv show`. If you don't specify a commit, it will search for the available commits. If you specify a commit it'll show you those results. For example:

```bash
$ asv show -b "views"
Commits with results:

Machine    : mimir.mobility.unimelb.net.au
Environment: conda-py3.7-h5py-memory_profiler-natsort-numpy-pandas-scipy

    61eb5bb7
    e9ccfc33
    22f12994
    0ebe187e
```

```bash
$ asv show -b "views" 0ebe187e
Commit: 0ebe187e <views-of-views>

views.SubsetMemorySuite.track_repeated_subset_memratio [mimir.mobility.unimelb.net.au/conda-py3.7-h5py-memory_profiler-natsort-numpy-pandas-scipy]
  ok
  ======= ======= ========== ============ ===================== ====================== ======================
  --                                                                   index_kind
  --------------------------------------- -------------------------------------------------------------------
   n_obs   n_var   attr_set   subset_dim         intarray             boolarray                slice
  ======= ======= ========== ============ ===================== ====================== ======================
    100     100     X-csr        obs               2.84           1.7916666666666667            0.5
    100     100     X-csr        var        2.5357142857142856    1.8695652173913044     0.5652173913043478
    100     100    X-dense       obs        3.1739130434782608    1.6538461538461537            0.6
...
```

You can compare two commits with `asv compare`

```bash
$ asv compare e9ccfc 0ebe187e
All benchmarks:

       before           after         ratio
     [e9ccfc33]       [0ebe187e]
     <master>         <views-of-views>
-            2.16  1.7916666666666667     0.83  views.SubsetMemorySuite.track_repeated_subset_memratio(100, 100, 'X-csr', 'obs', 'boolarray')
+ 2.533333333333333             2.84     1.12  views.SubsetMemorySuite.track_repeated_subset_memratio(100, 100, 'X-csr', 'obs', 'intarray')
- 1.1923076923076923              0.5     0.42  views.SubsetMemorySuite.track_repeated_subset_memratio(100, 100, 'X-csr', 'obs', 'slice')
  1.9615384615384615  1.8695652173913044     0.95  views.SubsetMemorySuite.track_repeated_subset_memratio(100, 100, 'X-csr', 'var', 'boolarray')
```

### View in the browser:

You can view the benchmarks in the browser with `asv publish` followed by `asv preview`. If you want to include benchmarks of a local branch, I think you'll have to add that branch to the `"branches"` list in `asv.conf.json`.
