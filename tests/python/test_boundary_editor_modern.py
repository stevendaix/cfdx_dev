"""Modern Fluent-style boundary workflow; no solver capability is invented."""

import math

import pytest

from cfdx.boundary_editor import (
    BOUNDARY_PRESETS, apply_boundary_condition, boundary_runtime_support,
    condition_for_field,
)
from cfdx.boundary_setup import validate_boundary_definition
from cfdx.case import Case
from cfdx.validation import validate_case


def test_single_field_remains_legacy_wire_shape():
    result = apply_boundary_condition(
        {}, patch_type="inlet", field="pressure", mode="FIXED_VALUE", value=101325.0,
    )
    assert result == {
        "type": "inlet", "fields": ["pressure"],
        "scalar_type": "FIXED_VALUE", "value": 101325.0, "gradient": 0.0,
    }
    validate_boundary_definition(result)
    assert boundary_runtime_support(result)[0]


def test_switching_fields_stages_without_losing_the_old_scalar_condition():
    previous = apply_boundary_condition(
        {"type": "inlet", "vendor_extension": {"key": 7}},
        patch_type="inlet", field="pressure", mode="FIXED_VALUE", value=101325.0,
    )
    original = dict(previous)
    new = apply_boundary_condition(
        previous, patch_type="inlet", field="velocity", mode="FIXED_VALUE",
        velocity=(4.0, -1.0, 0.0),
    )
    assert previous == original
    assert new["vendor_extension"] == {"key": 7}
    assert new["fields"] == ["pressure", "velocity"]
    assert condition_for_field(new, "pressure")["value"] == 101325.0
    assert condition_for_field(new, "velocity")["velocity_value"] == [4.0, -1.0, 0.0]
    assert "value" not in new and "velocity_value" not in new
    assert not boundary_runtime_support(new)[0]
    case = Case(name="test")
    case.set_boundary("inlet", **new)
    errors = validate_case(case).errors
    assert any(d.code == "BOUNDARY_SCHEMA" and "multi-field" in d.message for d in errors)


def test_editing_staged_field_updates_only_that_field():
    old = apply_boundary_condition(
        {}, patch_type="outlet", field="pressure", mode="FIXED_VALUE", value=0.0
    )
    staged = apply_boundary_condition(
        old, patch_type="outlet", field="temperature", mode="FIXED_VALUE", value=300.0,
    )
    updated = apply_boundary_condition(
        staged, patch_type="outlet", field="temperature", mode="FIXED_VALUE", value=310.0,
    )
    assert condition_for_field(updated, "pressure")["value"] == 0.0
    assert condition_for_field(updated, "temperature")["value"] == 310.0
    assert condition_for_field(staged, "temperature")["value"] == 300.0


def test_type_change_requires_acknowledgement_and_clears_stale_conditions():
    original = apply_boundary_condition(
        {"vendor_extension": {"keep": True}},
        patch_type="wall", field="velocity", mode="FIXED_VALUE", velocity=(0, 0, 0),
    )
    with pytest.raises(ValueError, match="explicit confirmation"):
        apply_boundary_condition(original, patch_type="symmetry")
    updated = apply_boundary_condition(
        original, patch_type="symmetry", allow_type_change=True
    )
    assert updated == {"type": "symmetry", "vendor_extension": {"keep": True}}
    validate_boundary_definition(updated)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf"), True, "4.0"])
def test_invalid_values_rejected_without_mutating_previous(bad):
    previous = {"type": "inlet", "fields": ["pressure"],
                "scalar_type": "FIXED_VALUE", "value": 1.0, "gradient": 0.0}
    with pytest.raises(ValueError, match="finite"):
        apply_boundary_condition(
            previous, patch_type="inlet", field="pressure", mode="FIXED_VALUE", value=bad
        )
    assert previous["value"] == 1.0


def test_invalid_velocity_shape_and_mode_rejected():
    with pytest.raises(ValueError, match="XYZ"):
        apply_boundary_condition({}, patch_type="inlet", field="velocity", mode="FIXED_VALUE", velocity=[1, 2])
    with pytest.raises(ValueError, match="unsupported"):
        apply_boundary_condition({}, patch_type="inlet", field="velocity", mode="ROBIN")


def test_presets_never_claim_solver_commit():
    names = {p.name for p in BOUNDARY_PRESETS}
    assert {"Velocity inlet", "Pressure outlet", "No-slip wall", "Symmetry"} <= names
    assert all(p.patch_type and p.explanation for p in BOUNDARY_PRESETS)
