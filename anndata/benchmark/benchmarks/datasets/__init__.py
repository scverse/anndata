"""This module should hold the code to read in datasets I'll use for benchmarking."""
from pathlib import Path as _Path

DATASETS_PATH = _Path(__file__).parent.absolute()

from .immune_cell_atlas import ICA_BoneMarrow_full, ICA_BoneMarrow_Donor1
from .pbmc_10x_v3 import PBMC_10X_NextGEM

DATASETS = [ICA_BoneMarrow_full, ICA_BoneMarrow_Donor1, PBMC_10X_NextGEM]


def list_available(datasets=DATASETS):
    """List datasets currently setup."""
    return list(filter(lambda x: x.is_available(), datasets))
