from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from anndata.tests.helpers import subset_func  # noqa: F401

if TYPE_CHECKING:
    from types import EllipsisType


@pytest.fixture
def backing_h5ad(tmp_path):
    return tmp_path / "test.h5ad"


@pytest.fixture(
    params=[
        pytest.param((..., (slice(None), slice(None))), id="ellipsis"),
        pytest.param(((...,), (slice(None), slice(None))), id="ellipsis_tuple"),
        pytest.param(
            ((..., slice(0, 10)), (slice(None), slice(0, 10))), id="obs-ellipsis"
        ),
        pytest.param(
            ((slice(0, 10), ...), (slice(0, 10), slice(None))), id="var-ellipsis"
        ),
        pytest.param(
            ((slice(0, 10), slice(0, 10), ...), (slice(0, 10), slice(0, 10))),
            id="obs-var-ellipsis",
        ),
        pytest.param(
            ((..., slice(0, 10), slice(0, 10)), (slice(0, 10), slice(0, 10))),
            id="ellipsis-obs-var",
        ),
        pytest.param(
            ((slice(0, 10), ..., slice(0, 10)), (slice(0, 10), slice(0, 10))),
            id="obs-ellipsis-var",
        ),
    ]
)
def ellipsis_index_with_equivalent(
    request,
) -> tuple[tuple[EllipsisType | slice, ...] | EllipsisType, tuple[slice, slice]]:
    return request.param


@pytest.fixture
def ellipsis_index(
    ellipsis_index_with_equivalent: tuple[
        tuple[EllipsisType | slice, ...] | EllipsisType, tuple[slice, slice]
    ],
) -> tuple[EllipsisType | slice, ...] | EllipsisType:
    return ellipsis_index_with_equivalent[0]


@pytest.fixture
def equivalent_ellipsis_index(
    ellipsis_index_with_equivalent: tuple[
        tuple[EllipsisType | slice, ...] | EllipsisType, tuple[slice, slice]
    ],
) -> tuple[slice, slice]:
    return ellipsis_index_with_equivalent[1]
