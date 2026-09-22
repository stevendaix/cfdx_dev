import pytest
from cfdx.initialization import InitializationMode, InitializationSpec, supported_initialization_modes

def test_supported_initialization_modes_are_explicit():
    assert supported_initialization_modes()==(InitializationMode.UNIFORM,InitializationMode.FIELD)

def test_uniform_initialization_requires_field_and_value():
    InitializationSpec(InitializationMode.UNIFORM,"U",0.0).validate()
    with pytest.raises(ValueError): InitializationSpec(InitializationMode.UNIFORM,"U").validate()

def test_field_initialization_requires_source():
    InitializationSpec(InitializationMode.FIELD,"restart/U").validate()
    with pytest.raises(ValueError): InitializationSpec(InitializationMode.FIELD).validate()

def test_initialization_requires_restart():
    from cfdx.setup_model import ChangeImpact
    assert InitializationSpec(InitializationMode.UNIFORM,"U",1.0).change_impact is ChangeImpact.REQUIRES_RESTART
