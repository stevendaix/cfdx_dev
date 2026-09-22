import pytest
from cfdx.materials import MaterialSpec

def test_material_properties_have_physical_units():
    m=MaterialSpec("air",1.2,1.8e-5,1005,0.026)
    m.validate()
    units={p.name:p.unit for p in m.parameters()}
    assert units["density"]=="kg/m^3"
    assert units["conductivity"]=="W/(m K)"

def test_material_rejects_nonphysical_values():
    with pytest.raises(ValueError,match="density"): MaterialSpec("bad",0,1e-3).validate()
    with pytest.raises(ValueError,match="specific heat"): MaterialSpec("bad",1,1e-3,0).validate()