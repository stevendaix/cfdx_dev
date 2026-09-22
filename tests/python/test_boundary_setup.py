import pytest
from cfdx.boundary_setup import fields_for_field, validate_boundary_definition

def test_velocity_boundary_has_typed_components():
    fields={f.name:f for f in fields_for_field("velocity")}
    assert fields["value_x"].unit=="m/s"

def test_unknown_boundary_type_is_rejected():
    with pytest.raises(KeyError): validate_boundary_definition({"type":"not-supported"})

def test_velocity_boundary_rejects_unknown_mode():
    with pytest.raises(ValueError): validate_boundary_definition({"type":"inlet","velocity_type":"ROBIN"},known_fields=("velocity",))
