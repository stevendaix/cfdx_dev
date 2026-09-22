import pytest
from cfdx.boundary_setup import SCALAR_TYPES, boundary_defaults, scalar_fields_for_boundary_type

def test_scalar_boundary_contract_matches_cpp_types():
    assert SCALAR_TYPES==("FIXED_VALUE","ZERO_GRADIENT","FIXED_GRADIENT")
    fields={f.name:f for f in fields_for_boundary_type("inlet")}
    assert fields["value"].default==0.0
    assert fields["gradient"].kind.value=="real"

def test_boundary_defaults_are_type_specific():
    assert boundary_defaults("periodic")=={}
    assert set(boundary_defaults("wall"))=={"type","value","gradient"}
    with pytest.raises(KeyError): fields_for_boundary_type("bogus")