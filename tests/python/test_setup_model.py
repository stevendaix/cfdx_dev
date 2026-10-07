from cfdx.setup_model import ChangeImpact, Parameter, ParameterType, typed_parameter

def test_parameter_carries_change_impact():
    p=typed_parameter("mesh",2,impact=ChangeImpact.REQUIRES_REBUILD)
    assert p.impact is ChangeImpact.REQUIRES_REBUILD

def test_parameter_rejects_nonfinite_real():
    import math, pytest
    with pytest.raises(ValueError): Parameter("x",math.inf,ParameterType.REAL).validate()
