import pytest


@pytest.mark.gpu
def test_gpu():
    """
    For testing that the gpu mark works
    """
    import cupy  # This test shouldn't run if cupy isn't installed

    cupy.ones(1)
