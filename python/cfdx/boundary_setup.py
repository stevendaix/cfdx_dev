"""Patch-aware boundary field contracts mirroring supported C++ solver inputs."""
from __future__ import annotations
from dataclasses import dataclass
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
    fields: tuple[BoundaryFieldSpec,...]

SCALAR_TYPES=("FIXED_VALUE","ZERO_GRADIENT","FIXED_GRADIENT")
VELOCITY_TYPES=("FIXED_VALUE","ZERO_GRADIENT")

VECTOR_FIELDS=(
    BoundaryFieldSpec("type",ParameterType.CHOICE,None,"ZERO_GRADIENT"),
    BoundaryFieldSpec("value_x",ParameterType.REAL,"m/s",0.0),
    BoundaryFieldSpec("value_y",ParameterType.REAL,"m/s",0.0),
    BoundaryFieldSpec("value_z",ParameterType.REAL,"m/s",0.0),
)
SCALAR_FIELDS=(
    BoundaryFieldSpec("type",ParameterType.CHOICE,None,"ZERO_GRADIENT"),
    BoundaryFieldSpec("value",ParameterType.REAL,None,0.0),
    BoundaryFieldSpec("gradient",ParameterType.REAL,None,0.0),
)

BOUNDARY_FIELD_SPECS={
 "inlet": VECTOR_FIELDS,
 "outlet": VECTOR_FIELDS,
 "wall": VECTOR_FIELDS,
 "symmetry": (BoundaryFieldSpec("type",ParameterType.CHOICE,None,"ZERO_GRADIENT"),),
 "periodic": (),
 "interface": (),
 "empty": (),
}

FIELD_SPECS={
 "velocity": VECTOR_FIELDS,
 "pressure": SCALAR_FIELDS,
 "temperature": SCALAR_FIELDS,
 "k": SCALAR_FIELDS,
 "epsilon": SCALAR_FIELDS,
 "omega": SCALAR_FIELDS,
}

def fields_for_boundary_type(boundary_type: str):
    try: return BOUNDARY_FIELD_SPECS[boundary_type]
    except KeyError as exc: raise KeyError(f"unsupported boundary type: {boundary_type}") from exc

def fields_for_field(field: str):
    try: return FIELD_SPECS[field]
    except KeyError as exc: raise KeyError(f"unsupported boundary field: {field}") from exc

def validate_boundary_definition(boundary: dict, *, known_fields: tuple[str,...]=()) -> None:
    if not isinstance(boundary,dict): raise ValueError("boundary definition must be a mapping")
    kind=boundary.get("type")
    fields_for_boundary_type(kind)
    for field in known_fields:
        specs=fields_for_field(field)
        for spec in specs:
            if spec.name=="type":
                continue
            if spec.name in boundary and spec.kind is ParameterType.REAL and not isinstance(boundary[spec.name],(int,float)):
                raise ValueError(f"{field}.{spec.name} must be numeric")
    if kind in {"inlet","outlet","wall"} and "velocity" in known_fields:
        if boundary.get("velocity_type","ZERO_GRADIENT") not in VELOCITY_TYPES:
            raise ValueError("unsupported velocity boundary type")

def boundary_defaults(boundary_type: str) -> dict[str,object]:
    return {f.name:f.default for f in fields_for_boundary_type(boundary_type)}
