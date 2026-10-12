from cfdx.case import Case
from cfdx.setup_model import ChangeImpact, ParameterType
from cfdx.setup_schema import build_setup_schema, setup_fields_for_selection

def test_setup_schema_exposes_representative_domains_from_one_contract() -> None:
    case = Case()
    case.enable("incompressible")
    case.enable("energy")
    case.set_numerics(cfl=0.8)
    case.set_boundary("inlet", type="inlet")
    schema = build_setup_schema(case)
    assert [domain.key for domain in schema.domains] == ["physics", "boundaries", "numerics", "solver"]
    assert schema.field("numerics.cfl").kind is ParameterType.REAL
    assert schema.field("numerics.cfl").value(case) == 0.8
    assert schema.field("execution.solver").impact is ChangeImpact.REQUIRES_RESTART
    assert schema.domain("physics").fields

def test_setup_schema_reuses_existing_boundary_and_physics_capabilities() -> None:
    case = Case()
    case.enable("incompressible")
    case.set_boundary("inlet", type="inlet")
    schema = build_setup_schema(case)
    physics_keys = {field.key for field in schema.domain("physics").fields}
    boundary_keys = {field.key for field in schema.domain("boundaries").fields}
    assert "physics.incompressible.density" in physics_keys
    assert "boundaries.inlet.value_x" in boundary_keys
    assert "boundaries.inlet.value_y" in boundary_keys

def test_setup_fields_for_selection_is_headless_and_value_based() -> None:
    case = Case()
    case.set_numerics(cfl=1.2)
    fields = setup_fields_for_selection(case, "numerics")
    assert [field.key for field in fields] == ["numerics.cfl"]
    assert fields[0].value(case) == 1.2
