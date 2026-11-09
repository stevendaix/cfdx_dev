"""Field-aware patch matrix for the Workbench editor (no Qt dependencies).

The existing single-field Python editor contract is not yet a certified mapping
onto the C++ CaseSetup HDF5 attribute. Multi-field configurations are blocked.
Additional field conditions can be staged losslessly but must be visibly marked
as blocked; no silent projection to Fluent/OpenFOAM-specific BC semantics.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import isfinite

from .boundary_editor import boundary_runtime_support, condition_for_field
from .boundary_setup import FIELD_SPECS

FIELD_UNITS = {
    "velocity": "m/s",
    "pressure": "Pa",
    "temperature": "K",
    "k": "m^2/s^2",
    "epsilon": "m^2/s^3",
    "omega": "1/s",
}
_GEOMETRIC = {"symmetry", "periodic", "interface", "empty"}
_LEGACY = frozenset(
    ("velocity_type", "velocity_value", "scalar_type", "value", "gradient")
)


@dataclass(frozen=True, slots=True)
class PatchFieldRow:
    patch: str
    patch_type: str
    field: str
    mode: str
    value: str
    units: str
    native_supported: bool
    runtime_note: str


def field_rows(patch: str, boundary: Mapping[str, object]) -> tuple[PatchFieldRow, ...]:
    """Return one row per explicitly configured field, never inferred values."""
    if not isinstance(patch, str) or not patch:
        raise ValueError("patch name must be provided")
    kind = boundary.get("type")
    if not isinstance(kind, str):
        raise ValueError("invalid patch type")  # noqa: TRY004
    fields = boundary.get("fields", [])
    if (
        not isinstance(fields, (list, tuple))
        or any(not isinstance(name, str) or name not in FIELD_SPECS for name in fields)
        or len(fields) != len(set(fields))
    ):
        raise ValueError("invalid or duplicate boundary field names")
    allowed, reason = boundary_runtime_support(boundary)
    rows = []
    for field in fields:
        data = condition_for_field(boundary, field)
        mode = data.get("velocity_type" if field == "velocity" else "scalar_type")
        if not isinstance(mode, str):
            raise ValueError(f"missing {field} boundary condition type")  # noqa: TRY004
        if field == "velocity":
            vec = data.get("velocity_value")
            if not isinstance(vec, (list, tuple)) or len(vec) != 3:
                raise ValueError("velocity boundary must contain XYZ")
            if any(
                not isinstance(value, (int, float))
                or isinstance(value, bool)
                or not isfinite(value)
                for value in vec
            ):
                raise ValueError("non-finite velocity component")
            value = "(" + ", ".join(f"{float(v):.7g}" for v in vec) + ")"
        elif mode == "FIXED_GRADIENT":
            v = data.get("gradient")
            if (
                not isinstance(v, (int, float))
                or isinstance(v, bool)
                or not isfinite(v)
            ):
                raise ValueError("invalid boundary gradient")
            value = f"grad {float(v):.7g}"
        elif mode == "ZERO_GRADIENT":
            value = "∂/∂n = 0"
        else:
            v = data.get("value")
            if (
                not isinstance(v, (int, float))
                or isinstance(v, bool)
                or not isfinite(v)
            ):
                raise ValueError("invalid boundary value")
            value = f"{float(v):.7g}"
        rows.append(
            PatchFieldRow(
                patch, kind, field, mode, value, FIELD_UNITS[field], allowed, reason
            )
        )
    return tuple(rows)


def remove_field_condition(
    boundary: Mapping[str, object],
    field: str,
) -> dict[str, object]:
    """Delete one explicit field, collapse to the certified legacy wire shape.

    No mutation occurs before the new definition has been validated.
    Other patch metadata and remaining field values are preserved.
    """
    fields = boundary.get("fields", [])
    if (
        not isinstance(fields, list)
        or field not in fields
        or len(fields) != len(set(fields))
    ):
        raise ValueError("field is not configured on this patch")
    remaining = [name for name in fields if name != field]
    values = {name: condition_for_field(boundary, name) for name in remaining}
    if any(not values[name] for name in remaining):
        raise ValueError("cannot remove a field when another condition is incomplete")
    result = {
        key: value
        for key, value in boundary.items()
        if key not in _LEGACY | {"field_conditions", "fields"}
    }
    if remaining:
        result["fields"] = remaining
    if len(remaining) == 1:
        result.update(values[remaining[0]])
    elif len(remaining) > 1:
        result["field_conditions"] = values
    return result


def boundary_matrix(boundaries: Mapping[str, object]) -> tuple[PatchFieldRow, ...]:
    result = []
    for name in sorted(boundaries):
        data = boundaries[name]
        if not isinstance(data, Mapping):
            raise ValueError(f"invalid patch definition: {name}")  # noqa: TRY004
        result.extend(field_rows(name, data))
    return tuple(result)


def setup_advice(
    boundaries: Mapping[str, object], mesh_patch_names=()
) -> tuple[str, ...]:
    """Non-blocking engineering guidance, never substituted for a solver gate."""
    messages = []
    known = set(mesh_patch_names)
    if known:
        unknown = set(boundaries) - known
        missing = known - set(boundaries)
        if unknown:
            messages.append(
                "Configured patches absent from mesh: " + ", ".join(sorted(unknown))
            )
        if missing:
            messages.append(
                "Mesh patches without explicit editor settings: "
                + ", ".join(sorted(missing))
            )
    pressure_refs = []
    for name, boundary in boundaries.items():
        if not isinstance(boundary, Mapping):
            messages.append(f"{name}: invalid boundary object")
            continue
        kind = boundary.get("type")
        if kind in _GEOMETRIC:
            continue
        fields = boundary.get("fields", [])
        if not fields:
            messages.append(f"{name}: no explicit field condition")
        if "pressure" in fields:
            cond = condition_for_field(boundary, "pressure")
            if cond.get("scalar_type") == "FIXED_VALUE":
                pressure_refs.append(name)
        ok, reason = boundary_runtime_support(boundary)
        if not ok:
            messages.append(f"{name}: native run blocked — {reason}")
    if boundaries and not pressure_refs:
        messages.append(
            "No explicit fixed pressure reference is recorded; verify the solver's pressure gauge strategy"
        )
    return tuple(messages)
