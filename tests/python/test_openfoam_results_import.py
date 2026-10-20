from __future__ import annotations

from pathlib import Path

import numpy as np

from cfdx.io.adapters.openfoam import OpenFOAMAdapter
from cfdx.io.interfaces import ConversionResult
from cfdx.io.schema import SourceInfo


def test_openfoam_adapter_imports_latest_scalar_and_vector_fields(
    tmp_path, monkeypatch
):
    case = tmp_path / "case"
    mesh = case / "constant" / "polyMesh"
    mesh.mkdir(parents=True)
    (mesh / "owner").write_text("owner fixture")

    adapter = OpenFOAMAdapter()
    monkeypatch.setattr(adapter, "_parse_label_list", lambda _: [0, 0, 1])

    def fake_read_fields(case_path, time, *, cell_count):
        assert Path(case_path) == case
        assert time == "latest"
        assert cell_count == 2
        return {
            "p": {
                "class": "volScalarField",
                "time": 10.0,
                "internal": np.array([101325.0, 101300.0]),
            },
            "U": {
                "class": "volVectorField",
                "time": 10.0,
                "internal": np.array([[1.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
            },
            "R": {
                "class": "volSymmTensorField",
                "time": 10.0,
                "internal": np.zeros((2, 6)),
            },
        }

    import cfdx.io.openfoam_case_reader as reader

    monkeypatch.setattr(reader, "read_openfoam_time_fields", fake_read_fields)
    result = ConversionResult(source=SourceInfo())

    assert adapter.import_results(str(case), result)
    assert result.field_data_source == "OpenFOAM time 10"
    assert result.available_fields == {"10": ["R", "U", "p"]}
    assert len(result.scalar_fields) == 1
    assert result.scalar_fields[0][0] == "p"
    np.testing.assert_array_equal(result.scalar_fields[0][1], [101325.0, 101300.0])
    assert len(result.vec_fields) == 1
    assert result.vec_fields[0][0] == "U"
    np.testing.assert_array_equal(
        result.vec_fields[0][1],
        [[1.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )
    assert any(
        finding.category == "results"
        and finding.severity.value == "unsupported_nonblocking"
        for finding in result.gap_report.findings
    )


def test_openfoam_adapter_reports_missing_case_for_results(tmp_path):
    adapter = OpenFOAMAdapter()
    result = ConversionResult(source=SourceInfo())

    assert not adapter.import_results(str(tmp_path / "missing"), result)
    assert result.gap_report.has_blocking()



def test_openfoam_polyhedral_mesh_builds_native_face_cell_csr():
    adapter = OpenFOAMAdapter()
    adapter._mesh_points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    faces = [
        [0, 1, 2],  # shared internal face
        [0, 3, 1],
        [1, 3, 2],
        [2, 3, 0],
        [0, 2, 4],
        [2, 1, 4],
        [1, 0, 4],
    ]
    owner = [0, 0, 0, 0, 1, 1, 1]
    neighbour = [1]
    patches = [{"name": "walls", "type": "wall", "n_faces": 6, "start_face": 1}]

    mesh = adapter._build_cells_from_faces(faces, owner, neighbour, patches)
    csr = mesh["cfdx_csr"]

    assert csr["n_points"] == 5
    assert csr["n_faces"] == 7
    assert csr["n_cells"] == 2
    np.testing.assert_array_equal(csr["owner"], owner)
    np.testing.assert_array_equal(csr["neighbour"], [1, -1, -1, -1, -1, -1, -1])
    np.testing.assert_array_equal(csr["cell_offsets"], [0, 4, 8])
    assert len(csr["cell_faces"]) == 8
    assert csr["boundary_patch_assignment"] == [
        {"name": "walls", "type": "wall", "face_ids": [1, 2, 3, 4, 5, 6]}
    ]
    # The surface triangles are for visualization only; solver topology is
    # represented by the native polyhedral CSR, not by these triangles.
    assert mesh["cells"]["triangle"].shape == (7, 3)



def test_native_openfoam_csr_is_written_to_cfdx_hdf5(tmp_path):
    import h5py

    from cfdx.io.hdf5_writer import write_case_cfdx_h5

    adapter = OpenFOAMAdapter()
    adapter._mesh_points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    faces = [
        [0, 1, 2],
        [0, 3, 1],
        [1, 3, 2],
        [2, 3, 0],
        [0, 2, 4],
        [2, 1, 4],
        [1, 0, 4],
    ]
    mesh = adapter._build_cells_from_faces(
        faces,
        [0, 0, 0, 0, 1, 1, 1],
        [1],
        [{"name": "walls", "type": "wall", "n_faces": 6, "start_face": 1}],
    )
    result = ConversionResult(
        source=SourceInfo(solver="OpenFOAM", format="openfoam"),
        mesh={"points": mesh["cfdx_csr"]["points"], **mesh},
        scalar_fields=[("p", np.array([101325.0, 101300.0]))],
        available_fields={"10": ["p"]},
        field_data_source="OpenFOAM time 10",
    )
    output = tmp_path / "case.cfdx.h5"

    write_case_cfdx_h5(result, output)

    with h5py.File(output, "r") as h5:
        assert int(h5.attrs["n_points"]) == 5
        assert int(h5.attrs["n_faces"]) == 7
        assert int(h5.attrs["n_cells"]) == 2
        np.testing.assert_array_equal(h5["owner"][:], [0, 0, 0, 0, 1, 1, 1])
        np.testing.assert_array_equal(h5["neighbour"][:], [1, -1, -1, -1, -1, -1, -1])
        np.testing.assert_array_equal(h5["cell_offsets"][:], [0, 4, 8])
        np.testing.assert_array_equal(h5["fields/scalar/p"][:], [101325.0, 101300.0])
        assert "boundary_patches" in h5.attrs
        assert "walls:1:6:0" in h5.attrs["boundary_patches"]
        assert h5.attrs["field_data_source"] == "OpenFOAM time 10"
        assert h5.attrs["available_fields_json"] == '{"10": ["p"]}'
