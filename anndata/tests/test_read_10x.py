from pathlib import Path

import h5py
import pytest
import anndata as ad
from anndata.tests.helpers import assert_equal
import shutil


ROOT = Path(__file__).parent / "data" / "10x"


@pytest.mark.parametrize(
    ["mtx_path", "h5_path"],
    [
        pytest.param(
            ROOT / "1.2.0" / "filtered_gene_bc_matrices" / "hg19_chr21",
            ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5",
        ),
        pytest.param(
            ROOT / "3.0.0" / "filtered_feature_bc_matrix",
            ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5",
        ),
    ],
)
@pytest.mark.parametrize("prefix", [None, "prefix_"])
def test_read_10x(tmp_path, mtx_path, h5_path, prefix):
    if prefix is not None:
        # Build files named "prefix_XXX.xxx" in a temporary directory.
        mtx_path_orig = mtx_path
        mtx_path = tmp_path / "filtered_gene_bc_matrices_prefix"
        mtx_path.mkdir()
        for item in mtx_path_orig.iterdir():
            if item.is_file():
                shutil.copyfile(item, mtx_path / f"{prefix}{item.name}")

    mtx = ad.read_10x_mtx(mtx_path, var_names="gene_symbols", prefix=prefix)
    h5 = ad.read_10x_h5(h5_path)

    # Drop genome column for comparing v3
    if "3.0.0" in str(h5_path):
        h5.var.drop(columns="genome", inplace=True)

    # Check equivalence
    assert_equal(mtx, h5)

    # Test that it can be written:
    from_mtx_pth = tmp_path / "from_mtx.h5ad"
    from_h5_pth = tmp_path / "from_h5.h5ad"

    mtx.write(from_mtx_pth)
    h5.write(from_h5_pth)

    assert_equal(ad.read_h5ad(from_mtx_pth), ad.read_h5ad(from_h5_pth))


def test_read_10x_h5_v1():
    spec_genome_v1 = ad.read_10x_h5(
        ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5",
        genome="hg19_chr21",
    )
    nospec_genome_v1 = ad.read_10x_h5(
        ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5"
    )
    assert_equal(spec_genome_v1, nospec_genome_v1)


def test_read_10x_h5_v2_multiple_genomes():
    genome1_v1 = ad.read_10x_h5(
        ROOT / "1.2.0" / "multiple_genomes.h5",
        genome="hg19_chr21",
    )
    genome2_v1 = ad.read_10x_h5(
        ROOT / "1.2.0" / "multiple_genomes.h5",
        genome="another_genome",
    )
    # the test data are such that X is the same shape for both "genomes",
    # but the values are different
    assert (genome1_v1.X != genome2_v1.X).sum() > 0, (
        "loading data from two different genomes in 10x v2 format. "
        "should be different, but is the same. "
    )


def test_read_10x_h5():
    spec_genome_v3 = ad.read_10x_h5(
        ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5",
        genome="GRCh38_chr21",
    )
    nospec_genome_v3 = ad.read_10x_h5(ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5")
    assert_equal(spec_genome_v3, nospec_genome_v3)


def test_error_10x_h5_legacy(tmp_path):
    onepth = ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5"
    twopth = tmp_path / "two_genomes.h5"
    with h5py.File(onepth, "r") as one, h5py.File(twopth, "w") as two:
        one.copy("hg19_chr21", two)
        one.copy("hg19_chr21", two, name="hg19_chr21_copy")
    with pytest.raises(ValueError):
        ad.read_10x_h5(twopth)
    ad.read_10x_h5(twopth, genome="hg19_chr21_copy")


def test_error_missing_genome():
    legacy_pth = ROOT / "1.2.0" / "filtered_gene_bc_matrices_h5.h5"
    v3_pth = ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5"
    with pytest.raises(ValueError, match=r".*hg19_chr21.*"):
        ad.read_10x_h5(legacy_pth, genome="not a genome")
    with pytest.raises(ValueError, match=r".*GRCh38_chr21.*"):
        ad.read_10x_h5(v3_pth, genome="not a genome")


def test_10x_h5_feature_types():
    # Tests that gex option doesn't, say, make the function return None
    h5_pth = ROOT / "3.0.0" / "filtered_feature_bc_matrix.h5"
    assert_equal(
        ad.read_10x_h5(h5_pth, feature_types="Gene Expression"),
        ad.read_10x_h5(h5_pth, feature_types=["Gene Expression"]),
    )
