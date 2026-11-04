from __future__ import annotations

from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from cfdx import dat_io
from cfdx.application import Application
from cfdx.case_inspector import inspect_cfdx_file, read_hdf5_dataset_sample
from cfdx.case_io import save_case
from cfdx.gui_3d import PyVistaQtView
from cfdx.session import CFDXSession


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
    sample = read_hdf5_dataset_sample(path, "/points", max_values=30)
    assert len(sample.values) == 30
    assert sample.total_values == 150_000_000
    assert set(sample.values) == {"0.0"}


def test_dataset_sample_is_bounded_and_rejects_unsafe_or_expensive_reads(tmp_path):
    path = tmp_path / "samples.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset("values", data=np.arange(1_000_000, dtype=np.float64))
        h5.create_dataset("strings", data=["one", "two"])
        h5.create_dataset(
            "oversized_chunk",
            shape=(1_000_000,),
            dtype=np.float64,
            chunks=(200_000,),
        )
        h5["soft_values"] = h5py.SoftLink("/values")
        h5["external_values"] = h5py.ExternalLink("external.h5", "/values")

    sample = read_hdf5_dataset_sample(path, "/values", max_values=12)
    assert len(sample.values) == 12
    assert sample.total_values == 1_000_000
    assert sample.values[0] == "0.0"
    assert sample.values[-1] == "11.0"

    with pytest.raises(ValueError, match="sampling dtype"):
        read_hdf5_dataset_sample(path, "/strings")
    with pytest.raises(ValueError, match="decompression sample budget"):
        read_hdf5_dataset_sample(path, "/oversized_chunk")
    with pytest.raises(ValueError, match="safe HDF5 path"):
        read_hdf5_dataset_sample(path, "/../outside")
    with pytest.raises(ValueError, match="between 1 and 256"):
        read_hdf5_dataset_sample(path, "/values", max_values=512)
    with pytest.raises(ValueError, match="internal HDF5 hard links"):
        read_hdf5_dataset_sample(path, "/soft_values")
    with pytest.raises(ValueError, match="internal HDF5 hard links"):
        read_hdf5_dataset_sample(path, "/external_values")


def test_inspect_case_written_by_canonical_case_io(tmp_path):
    path = tmp_path / "roundtrip.cfdx.h5"
    session = CFDXSession()

    save_case(session, path)
    summary = Application().inspect_case_file(path)

    assert summary.is_cfdx
    assert summary.has_project_config
    assert summary.case_name == session.case.name
    assert summary.physics == tuple(sorted(session.case.physics))
    assert summary.numerics == tuple(sorted(session.case.numerics))
    assert summary.boundaries == tuple(sorted(session.case.boundaries))
    assert not summary.has_mesh
    assert not summary.preview_supported


def test_inspect_converter_mesh_layout_without_case_config(tmp_path):
    path = tmp_path / "mesh.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_points"] = 4
        h5.attrs["n_faces"] = 4
        h5.attrs["n_cells"] = 1
        h5.attrs["boundary_patches"] = "wall:0:4:0"
        h5.attrs["source_case_name"] = "converted-channel"
        h5.attrs["case_setup_json"] = (
            '{"name":"converted-channel","physics":{"energy":true},'
            '"numerics":{"scheme":"SIMPLE"},"boundaries":{"wall":{}}}'
        )
        h5.create_dataset("points", data=np.zeros((4, 3)))
        h5.create_dataset("owner", data=[0, 0, 0, 0])
        h5.create_dataset("neighbour", data=[-1, -1, -1, -1])
        h5.create_dataset("face_vertices", data=[0, 1, 1, 2, 2, 3, 3, 0])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])

    summary = inspect_cfdx_file(path)

    assert summary.is_cfdx
    assert summary.has_mesh
    assert summary.has_case_config
    assert not summary.has_project_config
    assert summary.case_name == "converted-channel"
    assert summary.physics == ("energy",)
    assert summary.numerics == ("scheme",)
    assert summary.preview_supported
    assert summary.n_cells == 1

    application = Application()
    application.select("patch:stale", "mesh", "Stale patch")
    catalog = application.load_mesh(path)
    assert catalog.n_cells == 1
    assert application.selection.stable_id is None


def test_soft_linked_case_config_is_not_reported_as_project_config(tmp_path):
    path = tmp_path / "soft-linked-config.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX"
        h5.attrs["case_setup_json"] = (
            '{"name":"imported","physics":{"energy":true},'
            '"numerics":{"scheme":"SIMPLE"},"boundaries":{}}'
        )
        h5.create_dataset(
            "metadata/decoy_config",
            data='{"name":"not-canonical","physics":{},"numerics":{},"boundaries":{}}',
        )
        h5.create_group("case")
        h5["case/config"] = h5py.SoftLink("/metadata/decoy_config")

    summary = inspect_cfdx_file(path)

    assert summary.has_case_config
    assert not summary.has_project_config
    assert summary.case_name == "imported"


def test_renderer_rejects_stale_inspection_summary_before_loading_arrays(
    tmp_path,
):
    path = tmp_path / "mesh-mutated-after-inspection.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.create_dataset("points", shape=(4, 3), dtype="f8")
        h5.create_dataset("owner", shape=(4,), dtype="i8")

    summary = inspect_cfdx_file(path)
    assert summary.file_signature is not None

    with h5py.File(path, "a") as h5:
        h5.create_dataset("metadata/mutation", data=np.arange(4096, dtype=np.int64))

    from cfdx.gui_3d import PyVistaQtView

    with pytest.raises(ValueError, match="changed since metadata inspection"):
        PyVistaQtView.load_cfdx_mesh(None, str(path), summary=summary)
    with pytest.raises(ValueError, match="changed since metadata inspection"):
        PyVistaQtView.load_cfdx_dat(None, str(path), "unused.dat", summary=summary)



@pytest.mark.parametrize("missing_dataset", ["cell_faces", "cell_offsets"])
def test_dat_preview_reports_missing_cell_topology(
    tmp_path, monkeypatch, missing_dataset
):
    path = tmp_path / f"missing-{missing_dataset}.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.attrs["boundary_patches"] = "wall:0:4:0"
        h5.create_dataset("points", data=np.zeros((4, 3)))
        h5.create_dataset("owner", data=np.zeros(4, dtype=np.int64))
        h5.create_dataset("neighbour", data=-np.ones(4, dtype=np.int64))
        h5.create_dataset("face_vertices", data=[0, 1, 1, 2, 2, 3, 3, 0])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])
        h5.create_dataset("cell_faces", data=[0, 1, 2, 3])
        h5.create_dataset("cell_offsets", data=[0, 4])
        del h5[missing_dataset]

    summary = inspect_cfdx_file(path)
    assert summary.preview_supported

    monkeypatch.setattr(
        dat_io,
        "read_dat_restart",
        lambda _path: SimpleNamespace(cells=1, fields={}),
    )
    with pytest.raises(TypeError, match=f"requires an internal HDF5 {missing_dataset}"):
        PyVistaQtView.load_cfdx_dat(None, str(path), "unused.dat", summary=summary)



def test_dat_preview_rejects_oversized_cell_connectivity_before_reading(
    tmp_path, monkeypatch
):
    path = tmp_path / "oversized-cell-connectivity.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.attrs["boundary_patches"] = "wall:0:4:0"
        h5.create_dataset("points", data=np.zeros((4, 3)))
        h5.create_dataset("owner", data=np.zeros(4, dtype=np.int64))
        h5.create_dataset("neighbour", data=-np.ones(4, dtype=np.int64))
        h5.create_dataset("face_vertices", data=[0, 1, 1, 2, 2, 3, 3, 0])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])
        h5.create_dataset(
            "cell_faces",
            shape=(2_000_001,),
            dtype=np.int64,
            chunks=(4096,),
        )
        h5.create_dataset("cell_offsets", data=[0, 2_000_001])

    summary = inspect_cfdx_file(path)
    assert summary.preview_supported

    from cfdx import dat_io
    from cfdx.gui_3d import PyVistaQtView

    def unexpected_restart_read(_path):
        raise AssertionError(
            "DAT checkpoint must not be read after failed mesh preflight"
        )

    monkeypatch.setattr(dat_io, "read_dat_restart", unexpected_restart_read)
    with pytest.raises(ValueError, match="2M-entry eager-load budget"):
        PyVistaQtView.load_cfdx_dat(None, str(path), "unused.dat", summary=summary)


def test_inspect_does_not_preview_mesh_with_external_region_catalog(tmp_path):
    path = tmp_path / "external-regions.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.create_dataset("points", shape=(4, 3), dtype="f8")
        h5.create_dataset("owner", data=[0, 0, 0, 0])
        h5.create_dataset("neighbour", data=[-1, -1, -1, -1])
        h5.create_dataset("face_vertices", data=[0, 1, 1, 2, 2, 3, 3, 0])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])
        h5["regions"] = h5py.SoftLink("/missing")

    summary = inspect_cfdx_file(path)

    assert summary.has_mesh
    assert not summary.preview_supported
    assert "region catalog" in summary.preview_block_reason


def test_inspect_does_not_preview_when_region_group_contains_soft_links(tmp_path):
    path = tmp_path / "soft-region-entry.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.create_dataset("points", shape=(4, 3), dtype="f8")
        h5.create_dataset("owner", data=[0, 0, 0, 0])
        h5.create_dataset("neighbour", data=[-1, -1, -1, -1])
        h5.create_dataset("face_vertices", data=[0, 1, 1, 2, 2, 3, 3, 0])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])
        regions = h5.create_group("regions")
        regions["linked"] = h5py.SoftLink("/missing")

    summary = inspect_cfdx_file(path)

    assert summary.has_mesh
    assert summary.regions_truncated
    assert not summary.preview_supported
    assert "non-internal links" in summary.preview_block_reason


def test_inspect_refuses_small_mesh_with_huge_surface_connectivity(tmp_path):
    path = tmp_path / "oversized-connectivity.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.create_dataset("points", data=np.zeros((4, 3)))
        h5.create_dataset("owner", data=np.zeros(4, dtype=np.int64))
        h5.create_dataset("neighbour", data=-np.ones(4, dtype=np.int64))
        h5.create_dataset(
            "face_vertices",
            shape=(2_000_001,),
            dtype=np.int64,
            chunks=(4096,),
        )
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])

    summary = inspect_cfdx_file(path)

    assert not summary.preview_supported
    assert "connectivity" in summary.preview_block_reason

    from cfdx.gui_3d import PyVistaQtView

    with pytest.raises(ValueError, match="preview deferred"):
        PyVistaQtView.load_cfdx_mesh(None, str(path), summary=summary)
    with pytest.raises(ValueError, match="preview deferred"):
        PyVistaQtView.load_cfdx_dat(None, str(path), "unused.dat", summary=summary)


def test_preview_rejects_non_integer_surface_connectivity(tmp_path):
    path = tmp_path / "invalid-connectivity-dtype.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.create_dataset("points", data=np.zeros((4, 3)))
        h5.create_dataset("owner", data=np.zeros(4, dtype=np.float64))
        h5.create_dataset("neighbour", data=-np.ones(4, dtype=np.int64))
        h5.create_dataset("face_vertices", data=[0, 1, 1, 2, 2, 3, 3, 0])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])

    summary = inspect_cfdx_file(path)

    assert not summary.preview_supported
    assert "integer dtypes" in summary.preview_block_reason


def test_preview_rejects_externally_stored_topology(tmp_path):
    path = tmp_path / "external-topology.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = 1
        h5.create_dataset(
            "points",
            shape=(4, 3),
            dtype=np.float64,
            external=[("external-points.raw", 0, h5py.h5f.UNLIMITED)],
        )
        h5.create_dataset("owner", data=np.zeros(4, dtype=np.int64))
        h5.create_dataset("neighbour", data=-np.ones(4, dtype=np.int64))
        h5.create_dataset("face_vertices", data=[0, 1, 1, 2, 2, 3, 3, 0])
        h5.create_dataset("face_offsets", data=[0, 2, 4, 6, 8])
        h5.create_dataset("patch_face_ids", data=[0, 1, 2, 3])
        h5.create_dataset("patch_face_offsets", data=[0, 4])

    summary = inspect_cfdx_file(path)

    assert not summary.preview_supported
    assert "externally stored" in summary.preview_block_reason


def test_inspect_bounds_dataset_tree(tmp_path):
    path = tmp_path / "many.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        for index in range(10):
            h5.create_dataset(f"fields/scalar/f{index}", data=[float(index)])

    summary = inspect_cfdx_file(path, max_datasets=3)

    assert len(summary.datasets) == 3
    assert summary.datasets_truncated


def test_inspect_bounds_soft_link_scans_in_region_and_field_groups(tmp_path):
    path = tmp_path / "many-soft-links.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        regions = h5.create_group("regions")
        scalar_fields = h5.create_group("fields/scalar")
        for index in range(1_205):
            link_name = f"link_{index:04d}"
            regions[link_name] = h5py.SoftLink("/missing")
            scalar_fields[link_name] = h5py.SoftLink("/missing")

    summary = inspect_cfdx_file(path)

    assert summary.regions == ()
    assert summary.regions_truncated
    assert summary.scalar_fields == ()
    assert summary.scalar_fields_truncated


def test_inspect_bounds_hdf5_group_only_trees(tmp_path):
    path = tmp_path / "many-groups.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        for index in range(1_005):
            h5.create_group(f"metadata/group_{index:04d}")

    summary = inspect_cfdx_file(path, max_datasets=1)

    assert summary.is_cfdx
    assert summary.datasets == ()
    assert summary.tree_truncated


def test_inspect_bounds_root_attribute_payloads_and_count(tmp_path):
    oversized_path = tmp_path / "oversized-attribute.h5"
    with h5py.File(oversized_path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["oversized"] = "x" * 8192

    oversized_summary = inspect_cfdx_file(oversized_path)
    assert dict(oversized_summary.attributes)["oversized"].startswith("<omitted:")

    path = tmp_path / "many-attributes.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        for index in range(105):
            h5.attrs[f"metadata_{index:03d}"] = index

    summary = inspect_cfdx_file(path)

    assert summary.is_cfdx
    assert summary.attributes_truncated
    assert len(summary.attributes) == 100


def test_inspect_falls_back_to_cell_offsets_for_invalid_count_attribute(tmp_path):
    path = tmp_path / "invalid-count.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = "not-a-count"
        h5.create_dataset("points", shape=(8, 3), dtype="f8")
        h5.create_dataset("owner", shape=(6,), dtype="u8")
        h5.create_dataset("cell_offsets", data=[0, 2, 4, 6, 8])

    summary = inspect_cfdx_file(path)

    assert summary.has_mesh
    assert summary.n_cells == 4


def test_inspect_marks_region_totals_as_partial_when_catalog_is_truncated(tmp_path):
    path = tmp_path / "many-regions.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        for index in range(301):
            group = h5.create_group(f"regions/region_{index:03d}")
            group.attrs["n_cells"] = 1
            group.create_dataset("points", shape=(4, 3), dtype="f8")
            group.create_dataset("owner", shape=(4,), dtype="u8")
            group.create_dataset("cell_offsets", data=[0, 4])

    summary = inspect_cfdx_file(path)

    assert summary.has_mesh
    assert len(summary.regions) == 300
    assert summary.regions_truncated
    assert not summary.mesh_counts_complete
    assert summary.n_cells == 300


def test_inspect_bounds_field_catalogs(tmp_path):
    path = tmp_path / "many-fields.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        for index in range(305):
            h5.create_dataset(f"fields/scalar/p{index:03d}", shape=(1,), dtype="f8")

    summary = inspect_cfdx_file(path)

    assert len(summary.scalar_fields) == 300
    assert summary.scalar_fields_truncated
    with pytest.raises(ValueError, match="between 1 and 5000"):
        inspect_cfdx_file(path, max_datasets=10_000)


def test_inspect_rejects_non_hdf5_file(tmp_path):
    path = tmp_path / "not-a-case.cfdx.h5"
    path.write_text("not hdf5")

    with pytest.raises(OSError):
        inspect_cfdx_file(path)
