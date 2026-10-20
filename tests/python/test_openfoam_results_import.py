from __future__ import annotations

from pathlib import Path

import numpy as np

from cfdx.io.adapters.openfoam import OpenFOAMAdapter
from cfdx.io.interfaces import ConversionResult
from cfdx.io.schema import SourceInfo


def test_openfoam_adapter_imports_latest_scalar_and_vector_fields(tmp_path, monkeypatch):
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
        finding.category == "results" and finding.severity.value == "unsupported_nonblocking"
        for finding in result.gap_report.findings
    )


def test_openfoam_adapter_reports_missing_case_for_results(tmp_path):
    adapter = OpenFOAMAdapter()
    result = ConversionResult(source=SourceInfo())

    assert not adapter.import_results(str(tmp_path / "missing"), result)
    assert result.gap_report.has_blocking()
