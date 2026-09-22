from cfdx.physics_setup import TURBULENCE_MODELS, default_parameters, physics_spec
from cfdx.setup_model import ParameterType

def test_physics_schema_matches_supported_cpp_control_families():
    energy=physics_spec("energy")
    assert energy.requires==("incompressible",)
    fields={f.name:f for f in energy.fields}
    assert fields["density"].unit=="kg/m^3"
    assert fields["conductivity"].unit=="W/(m K)"

def test_turbulence_models_match_cpp_enum():
    params={p.name:p for p in default_parameters(physics_spec("turbulence"))}
    assert params["model"].kind is ParameterType.CHOICE
    assert params["model"].choices==TURBULENCE_MODELS
    assert params["molecular_viscosity"].unit=="Pa s"