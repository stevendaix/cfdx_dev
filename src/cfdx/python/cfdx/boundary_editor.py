"""Lossless, solver-aware boundary editor actions (no Qt dependency).

The current CFDX case contract has a legacy *single-field* encoding. Editing
another field must not silently erase the first one. Multiple independent
field conditions may be staged and saved in the case, but native execution is
blocked by validate_case until a multi-field solver mapping is certified.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import isfinite

from .boundary_setup import (
    BOUNDARY_FIELD_SPECS,
    FIELD_SPECS,
    SCALAR_TYPES,
    VELOCITY_TYPES,
)

_LEGACY_KEYS = frozenset(
    (
        "scalar_type",
        "value",
        "gradient",
        "velocity_type",
        "velocity_value",
    )
)
_NO_FIELDS = frozenset(("symmetry", "periodic", "interface", "empty"))


@dataclass(frozen=True, slots=True)
class BoundaryPreset:
    name: str
    patch_type: str
    field: str | None
    condition: str | None
    explanation: str


BOUNDARY_PRESETS = (
    BoundaryPreset(
        "Velocity inlet",
        "inlet",
        "velocity",
        "FIXED_VALUE",
        "Prescribe the velocity components; confirm magnitude and flow direction.",
    ),
    BoundaryPreset(
        "Pressure outlet",
        "outlet",
        "pressure",
        "FIXED_VALUE",
        "Prescribe an outlet pressure; verify the pressure reference and possible backflow.",
    ),
    BoundaryPreset(
        "No-slip wall",
        "wall",
        "velocity",
        "FIXED_VALUE",
        "Set wall velocity to zero for a stationary no-slip wall.",
    ),
    BoundaryPreset(
        "Moving wall",
        "wall",
        "velocity",
        "FIXED_VALUE",
        "Set a physically justified, nonzero wall velocity.",
    ),
    BoundaryPreset(
        "Symmetry",
        "symmetry",
        None,
        None,
        "Geometric symmetry; do not invent scalar or velocity values.",
    ),
)


def _number(value: object, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not isfinite(value)
    ):
        raise ValueError(f"{name} must be a finite real number")
    return float(value)


def _mode(field: str, mode: str) -> str:
    if field not in FIELD_SPECS:
        raise ValueError(f"unsupported field: {field!r}")
    allowed = (
        VELOCITY_TYPES
        if field == "velocity"
        else SCALAR_TYPES
        + (
            ("WALL_K_ZERO",)
            if field == "k"
            else ("WALL_OMEGA_SST",)
            if field == "omega"
            else ()
        )
    )
    if mode not in allowed:
        raise ValueError(f"unsupported {field} condition: {mode!r}")
    return mode


def condition_for_field(
    boundary: Mapping[str, object], field: str
) -> dict[str, object]:
    """Return stored condition for exactly one selected field, if present."""
    if field not in FIELD_SPECS:
        raise ValueError(f"unsupported field: {field!r}")
    staged = boundary.get("field_conditions", {})
    if isinstance(staged, dict) and field in staged and isinstance(staged[field], dict):
        return dict(staged[field])
    recorded = boundary.get("fields")
    if isinstance(recorded, list) and recorded == [field]:
        return {key: boundary[key] for key in _LEGACY_KEYS if key in boundary}
    return {}


def _canonical_condition(
    field: str,
    mode: str,
    *,
    value: object = 0.0,
    gradient: object = 0.0,
    velocity: object = (0.0, 0.0, 0.0),
) -> dict[str, object]:
    mode = _mode(field, mode)
    if field == "velocity":
        if not isinstance(velocity, (list, tuple)) or len(velocity) != 3:
            raise ValueError("velocity must be an XYZ vector")
        vector = [_number(v, "velocity component") for v in velocity]
        return {"velocity_type": mode, "velocity_value": vector}
    return {
        "scalar_type": mode,
        "value": _number(value, "boundary value"),
        "gradient": _number(gradient, "boundary gradient"),
    }


def apply_boundary_condition(
    existing: Mapping[str, object] | None,
    *,
    patch_type: str,
    field: str | None = None,
    mode: str | None = None,
    value: object = 0.0,
    gradient: object = 0.0,
    velocity: object = (0.0, 0.0, 0.0),
    allow_type_change: bool = False,
) -> dict[str, object]:
    """Build a new case boundary without mutating the previous one.

    For a single field the legacy schema is maintained byte-for-byte in shape.
    A second field is staged losslessly; validate_case reports the unimplemented
    native multi-field mapping as a blocking error.
    """
    if patch_type not in BOUNDARY_FIELD_SPECS:
        raise ValueError(f"unsupported patch type: {patch_type!r}")
    previous = dict(existing or {})
    old_type = previous.get("type")
    if old_type is not None and old_type != patch_type and not allow_type_change:
        raise ValueError("changing the patch kind requires explicit confirmation")
    if old_type is not None and old_type != patch_type:
        for key in (*_LEGACY_KEYS, "fields", "field_conditions"):
            previous.pop(key, None)
    result = dict(previous)
    result["type"] = patch_type
    if patch_type in _NO_FIELDS:
        for key in (*_LEGACY_KEYS, "fields", "field_conditions"):
            result.pop(key, None)
        return result
    if field is None or mode is None:
        raise ValueError("a field and condition are required for this patch type")
    condition = _canonical_condition(
        field, mode, value=value, gradient=gradient, velocity=velocity
    )
    old_fields = previous.get("fields", [])
    if not isinstance(old_fields, list) or any(
        not isinstance(item, str) for item in old_fields
    ):
        raise ValueError("existing boundary fields are invalid")
    if len(old_fields) != len(set(old_fields)):
        raise ValueError("duplicate boundary fields")
    updated_fields = list(old_fields)
    if field not in updated_fields:
        updated_fields.append(field)
    if len(updated_fields) == 1:
        result.pop("field_conditions", None)
        result["fields"] = [field]
        for key in _LEGACY_KEYS:
            result.pop(key, None)
        result.update(condition)
        return result
    staged = previous.get("field_conditions", {})
    if not isinstance(staged, dict):
        raise ValueError("existing field conditions must be a mapping")  # noqa: TRY004
    conditions = {
        name: dict(entry) for name, entry in staged.items() if isinstance(entry, dict)
    }
    for name in old_fields:
        if name not in conditions:
            legacy = condition_for_field(previous, name)
            if not legacy:
                raise ValueError(f"missing previous condition for {name}")
            conditions[name] = legacy
    conditions[field] = condition
    for key in _LEGACY_KEYS:
        result.pop(key, None)
    result["fields"] = updated_fields
    result["field_conditions"] = conditions
    return result


def boundary_runtime_support(boundary: Mapping[str, object]) -> tuple[bool, str]:
    """Fail closed on staged multi-field conditions without a native mapping."""
    fields = boundary.get("fields", [])
    if not isinstance(fields, list):
        return False, "invalid field list"
    if len(fields) > 1 or boundary.get("field_conditions"):
        return (
            False,
            "multi-field boundaries require explicit native CaseSetup export before CLI execution",
        )
    return (
        True,
        "single-field editor contract (native HDF5 CaseSetup bridge not certified)",
    )
