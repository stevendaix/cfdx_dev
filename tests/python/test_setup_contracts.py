from dataclasses import dataclass
import pytest
from cfdx import Case
from cfdx.setup_model import MeshSelection, Parameter, ParameterType, typed_parameter
from cfdx.validation import validate_case

@dataclass(frozen=True)
class FakeMesh:
    n_cells: int
    patches: tuple[str, ...]

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
    with pytest.raises(ValueError): MeshSelection("patch", -1)
    with pytest.raises(ValueError): MeshSelection("triangle", 0)

def test_case_validation_reports_errors_and_warnings_without_qt():
    case = Case(name="pipe")
    case.numerics["cfl"] = -1.0
    case.set_boundary("inlet", type="inlet", value=1.0)
    report = validate_case(case, FakeMesh(10, ("wall",)))
    assert not report.ok
    assert {d.code for d in report.errors} == {"CFL", "ORPHAN_BOUNDARY"}
    assert {d.code for d in report.warnings} == {"SOLVER_UNSET"}

def test_valid_case_has_no_errors():
    case = Case(name="pipe")
    case.execution.solver = "cfdx_solver"
    case.numerics["cfl"] = 1.0
    case.set_boundary("inlet", type="inlet", value=1.0)
    report = validate_case(case, FakeMesh(10, ("inlet",)))
    assert report.ok
    assert not report.errors
