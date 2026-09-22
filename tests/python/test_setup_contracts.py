from dataclasses import dataclass

import pytest

from cfdx import Case
from cfdx.setup_model import MeshSelection, Parameter, ParameterType, typed_parameter
from cfdx.validation import validate_case


@dataclass(frozen=True)
class FakeMesh:
    n_cells: int
    patches: tuple[str, ...]


def _configured_case() -> Case:
    case = Case(name="pipe")
    case.execution.solver = "cfdx_solver"
    case.numerics["cfl"] = 1.0
    case.physics["incompressible"] = {"enabled": True}
    case.materials["water"] = {
        "density": 1000.0,
        "dynamic_viscosity": 0.001,
        "cp": 4180.0,
        "conductivity": 0.6,
    }
    return case


def test_typed_parameter_validates_type_unit_and_choices():
    p = Parameter("temperature", 300.0, ParameterType.REAL, "K")
    p.validate()
    assert p.unit == "K"
    p = typed_parameter("model", "k-epsilon", choices=("laminar", "k-epsilon"))
    assert p.kind is ParameterType.CHOICE
    with pytest.raises(ValueError, match="unsupported value"):
        Parameter("model", "bogus", ParameterType.CHOICE, choices=("laminar",)).validate()


def test_mesh_selection_is_stable_and_validated():
    assert MeshSelection("patch", 3, "inlet").name == "inlet"
    with pytest.raises(ValueError):
        MeshSelection("patch", -1)
    with pytest.raises(ValueError):
        MeshSelection("triangle", 0)


def test_case_validation_reports_errors_and_warnings_without_qt():
    case = _configured_case()
    case.numerics["cfl"] = -1.0
    case.set_boundary("inlet", type="inlet", value=1.0)
    report = validate_case(case, FakeMesh(10, ("wall",)))
    assert not report.ok
    assert {d.code for d in report.errors} == {"CFL", "ORPHAN_BOUNDARY"}
    assert not report.warnings


def test_valid_case_has_no_errors():
    case = _configured_case()
    case.set_boundary("inlet", type="inlet", fields=["pressure"], scalar_type="ZERO_GRADIENT", value=1.0, gradient=0.0)
    report = validate_case(case, FakeMesh(10, ("inlet",)))
    assert report.ok
    assert not report.errors


def test_case_validation_reports_missing_physics_and_materials():
    case = Case(name="pipe")
    case.execution.solver = "cfdx_solver"
    report = validate_case(case)
    assert any(d.code == "PHYSICS_UNSET" for d in report.errors)


def test_case_validation_reports_dependency_and_material_errors():
    case = Case(name="pipe")
    case.execution.solver = "cfdx_solver"
    case.physics["energy"] = {"enabled": True}
    case.materials["bad"] = {"density": -1.0, "dynamic_viscosity": 0.0, "cp": 0.0, "conductivity": 0.0}
    report = validate_case(case)
    assert any(d.code == "PHYSICS_DEPENDENCY" for d in report.errors)
    assert any(d.code == "MATERIAL_VALUE" for d in report.errors)


def test_case_validation_rejects_invalid_boundary_schema():
    case = _configured_case()
    case.boundaries["inlet"] = {
        "type": "inlet",
        "velocity_type": "ROBIN",
        "fields": ["velocity"],
    }
    report = validate_case(case)
    assert any(d.code == "BOUNDARY_SCHEMA" for d in report.errors)
