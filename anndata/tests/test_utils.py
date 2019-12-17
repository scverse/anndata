import pytest
import zarr
import h5py

from anndata._io.utils import report_read_key_on_error, AnnDataReadError


def test_key_error(tmp_path):
    @report_read_key_on_error
    def read_attr(group):
        raise NotImplementedError()

    z = zarr.group()
    z["X"] = [1, 2, 3]
    z.create_group("group")
    with pytest.raises(AnnDataReadError) as e:
        read_attr(z["X"])
    assert "'/X'" in str(e.value)
    with pytest.raises(AnnDataReadError) as e:
        read_attr(z["group"])
    assert "'/group'" in str(e.value)

    h5pth = tmp_path / "test.h5"
    with h5py.File(h5pth, mode="a") as f:
        f["X"] = [1, 2, 3]
        f.create_group("group")
        with pytest.raises(AnnDataReadError) as e:
            read_attr(f["X"])
        assert "'/X'" in str(e.value)
        with pytest.raises(AnnDataReadError) as e:
            read_attr(f["group"])
        assert "'/group'" in str(e.value)
