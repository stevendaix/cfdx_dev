from __future__ import annotations

import h5py
import numpy as np
import pytest

from cfdx.application import Application
from cfdx.case_inspector import inspect_cfdx_file


def test_inspect_canonical_case_and_mesh_without_loading_arrays(tmp_path):
    path = tmp_path / "channel.cfdx.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX"
        h5.attrs["schema_version"] = 1
        h5.attrs["n_cells"] = 4
        h5.create_dataset("points", shape=(50_000_000, 3), dtype="f8", chunks=(1024, 3))
        h5.create_dataset("owner", shape=(120_000_000,), dtype="u8", chunks=(4096,))
        h5.create_dataset("cell_offsets", data=np.arange(5, dtype=np.uint64))
        h5.create_dataset(
            "face_vertices", shape=(360_000_000,), dtype="u8", chunks=(4096,)
        )
        h5.create_dataset(
            "face_offsets", shape=(120_000_001,), dtype="u8", chunks=(4096,)
        )
        h5.create_dataset("patch_face_ids", shape=(10,), dtype="u8")
        h5.create_dataset("patch_face_offsets", data=[0, 10])
        h5.attrs["boundary_patches"] = "wall:0:10:0"
        h5.create_dataset(
            "case/config",
            data=(
                '{"name":"channel","physics":{"energy":true},'
                '"numerics":{"cfl":1.0},"boundaries":{"wall":{}}}'
            ),
        )
        h5.create_dataset("fields/scalar/p", shape=(4,), dtype="f8")
        h5.create_dataset("fields/vector/U", shape=(4, 3), dtype="f8")

    summary = Application().inspect_case_file(path)

    assert summary.is_cfdx
    assert summary.has_case_config
    assert summary.has_mesh
    assert not summary.preview_supported
    assert summary.n_points == 50_000_000
    assert summary.n_cells == 4
    assert summary.scalar_fields == ("p",)
    assert summary.vector_fields == ("U",)
    assert summary.case_name == "channel"
    assert summary.physics == ("energy",)
    assert summary.has_project_config
    assert summary.numerics == ("cfl",)
    assert summary.boundaries == ("wall",)
    assert any(
        item.path == "/points" and item.size_bytes == 1_200_000_000
        for item in summary.datasets
    )


def test_inspect_converter_mesh_layout_without_case_config(tmp_path):
    path = tmp_path / "mesh.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_points"] = 4
        h5.attrs["n_faces"] = 4
        h5.attrs["n_cells"] = 1
        h5.create_dataset("points", data=np.zeros((4, 3)))
        h5.create_dataset("owner", data=[0, 0, 0, 0])
        h5.create_dataset("neighbour", data=[-1, -1, -1, -1])
        h5.create_dataset("neighbour", data=[-1, -1, -1, -1])
        h5.create_dataset("face_vertices", data=[0, 1, 2, 3])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])

    summary = inspect_cfdx_file(path)

    assert summary.is_cfdx
    assert summary.has_mesh
    assert not summary.has_case_config
    assert summary.preview_supported
    assert summary.n_cells == 1


def test_inspect_bounds_dataset_tree(tmp_path):
    path = tmp_path / "many.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        for index in range(10):
            h5.create_dataset(f"fields/scalar/f{index}", data=[float(index)])

    summary = inspect_cfdx_file(path, max_datasets=3)

    assert len(summary.datasets) == 3
    assert summary.datasets_truncated


def test_inspect_rejects_non_hdf5_file(tmp_path):
    path = tmp_path / "not-a-case.cfdx.h5"
    path.write_text("not hdf5")

    with pytest.raises(OSError):
        inspect_cfdx_file(path)
