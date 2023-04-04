import re
import sys
import scanpy as sc

from pathlib import Path

sys.path.append(str(Path(__file__).parent))


def download_dataset():
    adata = sc.read_10x_mtx(
        # the directory with the `.mtx` file
        "filtered_gene_bc_matrices/hg19/",
        # use gene symbols for the variable names (variables-axis index)
        var_names="gene_symbols",
        # write a cache file for faster subsequent reading
        cache=True,
    )

    adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
    adata.write("pbmc3k_raw.h5ad", compression="gzip")


def main():
    from argparse import ArgumentParser

    # parser = ArgumentParser(description="Setup datasets for benchmarking.")
    # args = parser.parse_args()
    # pat = re.compile(args.pattern)

    download_dataset()


if __name__ == "__main__":
    main()
