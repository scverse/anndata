# conda deactivate

mamba env remove -yn anndata-min-deps-test
mamba create -yn anndata-min-deps-test "python=3.9"

PACKAGES=`python3 ci/scripts/min-deps.py pyproject.toml --extra dev test`

# conda activate anndata-min-deps-test
# conda run -n anndata-min-deps-test pip install cupy-cuda12x


echo Installing $PACKAGES
conda run -n anndata-min-deps-test pip install $PACKAGES
conda run -n anndata-min-deps-test pip install pytest-xdist # cupy-cuda12x
conda run -n anndata-min-deps-test pip install -e . --no-deps
conda run -n anndata-min-deps-test pytest -n auto

conda list -n anndata-min-deps-tests
