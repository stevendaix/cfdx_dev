import pytest

from cfdx.boundary_editor import apply_boundary_condition, boundary_runtime_support
from cfdx.boundary_matrix import (
    FIELD_UNITS, boundary_matrix, field_rows, remove_field_condition, setup_advice,
)
from cfdx.boundary_setup import validate_boundary_definition


def _scalar(field, value):
    return apply_boundary_condition({}, patch_type="outlet", field=field,
                                    mode="FIXED_VALUE", value=value)


def test_pressure_and_temperature_rows_have_explicit_units():
    pressure = _scalar("pressure", 101325.0)
    staged = apply_boundary_condition(pressure, patch_type="outlet", field="temperature",
                                      mode="FIXED_VALUE", value=300.0)
    rows = field_rows("exit", staged)
    assert [(r.field, r.value, r.units) for r in rows] == [
        ("pressure", "101325", "Pa"), ("temperature", "300", "K"),
    ]
    assert all(not row.native_supported for row in rows)
    assert "native run blocked" in " ".join(setup_advice({"exit": staged}))
    assert FIELD_UNITS["omega"] == "1/s"


def test_field_removal_is_transactional_and_collapses_to_native_wire_shape():
    pressure = _scalar("pressure", 101325)
    pressure["private_extension"] = {"x": [1]}
    staged = apply_boundary_condition(pressure, patch_type="outlet", field="temperature",
                                      mode="FIXED_VALUE", value=310.0)
    original = dict(staged)
    result = remove_field_condition(staged, "temperature")
    assert staged == original
    assert result == pressure
    assert boundary_runtime_support(result)[0]
    validate_boundary_definition(result)


def test_removing_last_field_does_not_invent_a_default():
    result = remove_field_condition(_scalar("pressure", 0), "pressure")
    assert result == {"type": "outlet"}
    assert any("no explicit field" in message for message in setup_advice({"out": result}))


def test_geometry_and_mesh_patch_coverage_advice():
    boundaries = {"sym": {"type": "symmetry"}, "out": _scalar("pressure", 0)}
    assert boundary_matrix(boundaries)[0].patch == "out"
    advice = setup_advice(boundaries, mesh_patch_names=("wall", "out"))
    assert any("Configured patches absent" in item for item in advice)
    assert any("without explicit editor settings" in item for item in advice)


@pytest.mark.parametrize("name", ["not-configured", "temperature"])
def test_removing_nonexistent_field_is_rejected(name):
    with pytest.raises(ValueError, match="not configured"):
        remove_field_condition(_scalar("pressure", 1), name)


def test_bad_staged_vector_never_treated_as_native_mapping():
    staged = {"type": "inlet", "fields": ["velocity", "pressure"],
              "field_conditions": {"velocity": {"velocity_type": "FIXED_VALUE",
                                               "velocity_value": [1, 2]}}}
    with pytest.raises(ValueError, match="velocity"):
        field_rows("inlet", staged)


def test_direct_launch_and_restart_block_staged_boundary_before_process(tmp_path):
    import sys
    from cfdx import CFDXSession, ExecutionController, SolverRunner

    original = _scalar("pressure", 101325.0)
    staged = apply_boundary_condition(
        original, patch_type="outlet", field="temperature", mode="FIXED_VALUE", value=298.0,
    )
    session = CFDXSession()
    session.case.boundaries["exit"] = staged
    controller = ExecutionController(session, SolverRunner([sys.executable, "-c", "pass"]))
    with pytest.raises(ValueError, match="exit.*multi-field"):
        controller.start()
    assert session.state.value == "CREATED"
    assert controller.runner._thread is None

    dat = tmp_path / "checkpoint.dat"
    dat.write_bytes(b"test checkpoint")
    with pytest.raises(ValueError, match="exit.*multi-field"):
        controller.restart(dat)
    assert session.state.value == "CREATED"
    assert controller.runner._thread is None
