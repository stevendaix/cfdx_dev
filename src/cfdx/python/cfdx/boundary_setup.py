"""Patch-aware boundary field contracts mirroring supported C++ solver inputs."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from .setup_model import ChangeImpact, ParameterType


@dataclass(frozen=True)
class BoundaryFieldSpec:
    name: str
    kind: ParameterType
    unit: str | None = None
    default: object = None
    impact: ChangeImpact = ChangeImpact.REQUIRES_RESTART


@dataclass(frozen=True)
class BoundaryTypeSpec:
    key: str
    fields: tuple[BoundaryFieldSpec, ...]


SCALAR_TYPES = ("FIXED_VALUE", "ZERO_GRADIENT", "FIXED_GRADIENT")
VELOCITY_TYPES = ("FIXED_VALUE", "ZERO_GRADIENT")

VECTOR_FIELDS = (
    BoundaryFieldSpec("type", ParameterType.CHOICE, None, "ZERO_GRADIENT"),
    BoundaryFieldSpec("value_x", ParameterType.REAL, "m/s", 0.0),
    BoundaryFieldSpec("value_y", ParameterType.REAL, "m/s", 0.0),
    BoundaryFieldSpec("value_z", ParameterType.REAL, "m/s", 0.0),
)
SCALAR_FIELDS = (
    BoundaryFieldSpec("type", ParameterType.CHOICE, None, "ZERO_GRADIENT"),
    BoundaryFieldSpec("value", ParameterType.REAL, None, 0.0),
    BoundaryFieldSpec("gradient", ParameterType.REAL, None, 0.0),
)

BOUNDARY_FIELD_SPECS = {
    "inlet": VECTOR_FIELDS,
    "outlet": VECTOR_FIELDS,
    "wall": VECTOR_FIELDS,
    "symmetry": (
        BoundaryFieldSpec("type", ParameterType.CHOICE, None, "ZERO_GRADIENT"),
    ),
    "periodic": (),
    "interface": (),
    "empty": (),
}

FIELD_SPECS = {
    "velocity": VECTOR_FIELDS,
    "pressure": SCALAR_FIELDS,
    "temperature": SCALAR_FIELDS,
    "k": SCALAR_FIELDS,
    "epsilon": SCALAR_FIELDS,
    "omega": SCALAR_FIELDS,
}


def fields_for_boundary_type(boundary_type: str):
    try:
        return BOUNDARY_FIELD_SPECS[boundary_type]
    except KeyError as exc:
        raise KeyError(f"unsupported boundary type: {boundary_type}") from exc


def fields_for_field(field: str):
    try:
        return FIELD_SPECS[field]
    except KeyError as exc:
        raise KeyError(f"unsupported boundary field: {field}") from exc


def validate_boundary_definition(
    boundary: dict, *, known_fields: tuple[str, ...] = ()
) -> None:
    """Validate the supported single-field wire shape and staged field data."""
    if not isinstance(boundary, dict):
        raise ValueError("boundary definition must be a mapping")
    kind = boundary.get("type")
    fields_for_boundary_type(kind)
    fields = boundary.get("fields", ())
    if not isinstance(fields, (list, tuple)) or any(
        not isinstance(field, str) for field in fields
    ):
        raise ValueError("boundary fields must be a list of names")
    if len(set(fields)) != len(fields):
        raise ValueError("duplicate boundary fields")
    if kind in {"symmetry", "periodic", "interface", "empty"} and fields:
        raise ValueError(f"{kind} must not declare scalar/velocity field values")
    declared = tuple(dict.fromkeys((*fields, *known_fields)))
    for field in declared:
        fields_for_field(field)
    from .boundary_editor import boundary_runtime_support, condition_for_field

    for field in declared:
        if field not in fields:
            continue
        condition = condition_for_field(boundary, field)
        mode_key = "velocity_type" if field == "velocity" else "scalar_type"
        mode = condition.get(mode_key, "ZERO_GRADIENT")
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
            raise ValueError(f"unsupported {field} boundary type: {mode!r}")
        if field == "velocity":
            vector = condition.get("velocity_value", (0.0, 0.0, 0.0))
            if not isinstance(vector, (list, tuple)) or len(vector) != 3:
                raise ValueError("velocity_value must contain three components")
            for item in vector:
                if (
                    isinstance(item, bool)
                    or not isinstance(item, (int, float))
                    or not isfinite(item)
                ):
                    raise ValueError("velocity_value must contain finite numbers")
        else:
            for key in ("value", "gradient"):
                value = condition.get(key, 0.0)
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not isfinite(value)
                ):
                    raise ValueError(f"{field}.{key} must be a finite number")
    if kind in {"inlet", "outlet", "wall"} and "velocity" in declared:
        if (
            boundary.get("velocity_type", "ZERO_GRADIENT") not in VELOCITY_TYPES
            and len(fields) <= 1
        ):
            raise ValueError("unsupported velocity boundary type")
    runtime_ok, explanation = boundary_runtime_support(boundary)
    if not runtime_ok:
        raise ValueError(explanation)


def boundary_defaults(boundary_type: str) -> dict[str, object]:
    return {f.name: f.default for f in fields_for_boundary_type(boundary_type)}


def scalar_fields_for_boundary_type(boundary_type: str):
    """Compatibility accessor for boundary field specifications."""
    return fields_for_boundary_type(boundary_type)
