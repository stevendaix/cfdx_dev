from cfdx.setup_model import ChangeImpact, Parameter, ParameterType, typed_parameter

def test_parameter_carries_change_impact():
    p=typed_parameter("mesh",2,impact=ChangeImpact.REQUIRES_REBUILD)
    assert p.impact is ChangeImpact.REQUIRES_REBUILD

def test_parameter_rejects_nonfinite_real():
    import math, pytest
    with pytest.raises(ValueError): Parameter("x",math.inf,ParameterType.REAL).validate()


from cfdx.setup_model import SetupDiagnostic


def test_setup_diagnostic_supports_shared_application_contract():
    diagnostic = SetupDiagnostic(
        "warning", "CFL_HIGH", "CFL is high", "numerics.cfl",
        domain="execution", source="solver", affected_object="case-1",
        remediation="Reduce CFL before starting the run.",
    )
    assert diagnostic.domain == "execution"
    assert diagnostic.source == "solver"
    assert diagnostic.affected_object == "case-1"
    assert diagnostic.remediation is not None


def test_setup_diagnostic_rejects_invalid_severity():
    import pytest
    with pytest.raises(ValueError, match="severity"):
        SetupDiagnostic("fatal", "X", "bad")


def test_setup_diagnostic_serializes_shared_fields():
    diagnostic = SetupDiagnostic("error", "X", "bad", "case.x", domain="results", source="solver")
    assert diagnostic.to_dict()["domain"] == "results"
    assert diagnostic.to_dict()["source"] == "solver"
