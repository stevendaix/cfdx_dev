from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cfdx.io import openfoam_case_reader as reader


class FakeCondition(dict):
    def as_dict(self):
        return dict(self)


class FakeField:
    def __init__(self, path: Path, dimensions, internal, boundary):
        self.path = path
        self.dimensions = dimensions
        self.internal_field = internal
        self.boundary_field = boundary

    def __fspath__(self):
        return str(self.path)

    def __str__(self):
        return str(self.path)


class FakeTime:
    def __init__(self, path: Path, time: float, fields):
        self.path = path
        self.name = path.name
        self.time = time
        self.fields = fields

    def __iter__(self):
        return iter(self.fields.values())

    def __getitem__(self, name):
        return self.fields[name]


class FakeCase:
    def __init__(self, path: Path):
        self.path = path
        self.control_dict = {"application": "simpleFoam"}
        self.times = [
            FakeTime(
                path / "0",
                0.0,
                {
                    "U": FakeField(
                        path / "0" / "U",
                        (0, 1, -1, 0, 0, 0, 0),
                        np.array([[1.0, 0.0, 0.0]]),
                        {"inlet": FakeCondition(type="fixedValue", value=[1, 0, 0])},
                    )
                },
            )
        ]

    def __iter__(self):
        return iter(self.times)

    def __getitem__(self, key):
        return next(time for time in self.times if time.name == str(key))


@pytest.fixture
def case_dir(tmp_path, monkeypatch):
    (tmp_path / "system").mkdir()
    (tmp_path / "system" / "controlDict").write_text("placeholder")
    (tmp_path / "constant").mkdir()
    (tmp_path / "0").mkdir()
    monkeypatch.setattr(reader, "_foamlib", lambda: FakeCase)
    return tmp_path


def test_inspect_openfoam_case_reports_times_fields_and_boundary(case_dir):
    summary = reader.inspect_openfoam_case(case_dir)

    assert summary.application == "simpleFoam"
    assert [item.time for item in summary.times] == [0.0]
    assert summary.times[0].fields[0].name == "U"
    assert summary.times[0].fields[0].dimensions == (0, 1, -1, 0, 0, 0, 0)
    assert summary.times[0].fields[0].boundary_types == {"inlet": "fixedValue"}
    assert "system/controlDict" in summary.dictionaries


def test_read_openfoam_field_returns_internal_and_boundary_values(case_dir):
    result = reader.read_openfoam_field(case_dir, 0, "U")

    assert result["name"] == "U"
    assert result["time"] == 0.0
    assert result["dimensions"] == (0, 1, -1, 0, 0, 0, 0)
    np.testing.assert_array_equal(result["internal"], [[1.0, 0.0, 0.0]])
    assert result["boundary"]["inlet"]["type"] == "fixedValue"


def test_inspect_rejects_missing_case_directory(tmp_path):
    with pytest.raises(NotADirectoryError):
        reader.inspect_openfoam_case(tmp_path / "missing")


def test_reader_gives_actionable_error_without_foamlib(tmp_path, monkeypatch):
    (tmp_path / "case").mkdir()

    def missing_dependency():
        raise reader.OpenFOAMDependencyError("install foamlib")

    monkeypatch.setattr(reader, "_foamlib", missing_dependency)
    with pytest.raises(reader.OpenFOAMDependencyError, match="install foamlib"):
        reader.inspect_openfoam_case(tmp_path / "case")
