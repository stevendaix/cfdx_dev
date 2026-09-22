"""Patch-aware boundary condition schemas for the application layer."""
from __future__ import annotations

from dataclasses import dataclass

from .setup_model import ParameterType

@dataclass(frozen=True)
class BoundaryFieldSpec:
    name: str
    kind: ParameterType
    unit: str | None = None
    default: object = None

@dataclass(frozen=True)
class BoundaryTypeSpec:
    key: str
    fields: tuple[BoundaryFieldSpec, ...]

SCALAR_TYPES=("FIXED_VALUE","ZERO_GRADIENT","FIXED_GRADIENT")

# The scalar contract mirrors cfdx::physics::ScalarBoundaryCondition.
SCALAR_FIELDS=(
 BoundaryFieldSpec("type",ParameterType.CHOICE,None,"ZERO_GRADIENT"),
 BoundaryFieldSpec("value",ParameterType.REAL,None,0.0),
 BoundaryFieldSpec("gradient",ParameterType.REAL,None,0.0),
)

BOUNDARY_FIELD_SPECS={
 "inlet": SCALAR_FIELDS,
 "outlet": SCALAR_FIELDS,
 "wall": SCALAR_FIELDS,
 "symmetry": (BoundaryFieldSpec("type",ParameterType.CHOICE,None,"ZERO_GRADIENT"),),
 "periodic": (),
 "interface": (),
 "empty": (),
}

def scalar_fields_for_boundary_type(boundary_type: str) -> tuple[BoundaryFieldSpec,...]:
    try:
        return BOUNDARY_FIELD_SPECS[boundary_type]
    except KeyError as exc:
        raise KeyError(f"unsupported boundary type: {boundary_type}") from exc

def boundary_defaults(boundary_type: str) -> dict[str,object]:
    return {field.name:field.default for field in scalar_fields_for_boundary_type(boundary_type)}
