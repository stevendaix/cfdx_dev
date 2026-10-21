from __future__ import annotations

import numpy as np
import pytest

from cfdx.io import openfoam_case_reader as reader


def _write_scalar_field(path, *, time="0", value=1.25):
    path.write_text(
        f"""FoamFile
{{
    version 2.0;
    format ascii;
    class volScalarField;
    location "{time}";
    object p;
}}
dimensions [0 2 -2 0 0 0 0];
internalField uniform {value};
boundaryField
{{
    inlet
    {{
        type fixedValue;
        value uniform {value};
    }}
}}
"""
    )


def test_inspect_openfoam_case_reports_numeric_times_and_fields(tmp_path):
    system = tmp_path / "system"
    system.mkdir()
    (system / "controlDict").write_text(
        """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object controlDict;
}
application simpleFoam;
"""
    )
    for time in ("0", "2", "10"):
        directory = tmp_path / time
        directory.mkdir()
        _write_scalar_field(directory / "p", time=time)

    summary = reader.inspect_openfoam_case(tmp_path)

    assert summary.application == "simpleFoam"
    assert [item.name for item in summary.times] == ["0", "2", "10"]
    assert [item.time for item in summary.times] == [0.0, 2.0, 10.0]
    assert summary.times[0].fields[0].name == "p"
    assert summary.times[0].fields[0].dimensions == (0, 2, -2, 0, 0, 0, 0)
    assert summary.times[0].fields[0].boundary_types == {"inlet": "fixedValue"}
    assert "system/controlDict" in summary.dictionaries


def test_read_openfoam_field_returns_internal_and_boundary_values(tmp_path):
    time_dir = tmp_path / "0"
    time_dir.mkdir()
    _write_scalar_field(time_dir / "p", value=4.5)

    result = reader.read_openfoam_field(tmp_path, 0, "p")

    assert result["name"] == "p"
    assert result["time"] == 0.0
    assert result["dimensions"] == (0, 2, -2, 0, 0, 0, 0)
    np.testing.assert_array_equal(result["internal"], np.asarray(4.5))
    assert result["boundary"]["inlet"]["type"] == "fixedValue"


def test_read_openfoam_field_expands_uniform_value_when_requested(tmp_path):
    from cfdx.io.openfoam_case_reader import read_openfoam_time_fields

    time_dir = tmp_path / "5"
    time_dir.mkdir()
    _write_scalar_field(time_dir / "p", time="5", value=8.0)

    fields = read_openfoam_time_fields(tmp_path, "latest", cell_count=3)

    np.testing.assert_array_equal(fields["p"]["internal"], [8.0, 8.0, 8.0])


def test_inspect_rejects_missing_case_directory(tmp_path):
    with pytest.raises(NotADirectoryError):
        reader.inspect_openfoam_case(tmp_path / "missing")


def test_native_reader_rejects_binary_field(tmp_path):
    time_dir = tmp_path / "0"
    time_dir.mkdir()
    (time_dir / "p").write_text(
        """FoamFile
{
    version 2.0;
    format binary;
    class volScalarField;
    object p;
}
"""
    )

    with pytest.raises(ValueError, match="Binary OpenFOAM files"):
        reader.read_openfoam_field(tmp_path, "latest", "p")
