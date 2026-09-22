"""Structured setup schemas for supported CFDX physics controls.

This is a UI/application schema, not a second solver implementation. Names and
defaults mirror the existing C++ control structures; numerical evaluation stays
in the C++ physics layer.
"""
from __future__ import annotations

from dataclasses import dataclass

from .setup_model import Parameter, ParameterType

@dataclass(frozen=True)
class FieldSpec:
    name: str
    kind: ParameterType
    unit: str | None = None
    default: object = None

@dataclass(frozen=True)
class PhysicsSpec:
    key: str
    label: str
    fields: tuple[FieldSpec, ...] = ()
    requires: tuple[str, ...] = ()

INCOMPRESSIBLE = PhysicsSpec("incompressible","Incompressible flow",(
    FieldSpec("density",ParameterType.REAL,"kg/m^3",1.0),
    FieldSpec("kinematic_viscosity",ParameterType.REAL,"m^2/s",1.0e-3),
    FieldSpec("algorithm",ParameterType.CHOICE,None,"SIMPLE"),
))
ENERGY = PhysicsSpec("energy","Energy",(
    FieldSpec("density",ParameterType.REAL,"kg/m^3",1.0),
    FieldSpec("cp",ParameterType.REAL,"J/(kg K)",1000.0),
    FieldSpec("conductivity",ParameterType.REAL,"W/(m K)",1.0),
    FieldSpec("dt",ParameterType.REAL,"s",0.0),
),("incompressible",))
TURBULENCE = PhysicsSpec("turbulence","Turbulence",(
    FieldSpec("model",ParameterType.CHOICE,None,"LAMINAR"),
    FieldSpec("density",ParameterType.REAL,"kg/m^3",1.0),
    FieldSpec("molecular_viscosity",ParameterType.REAL,"Pa s",1.0e-3),
    FieldSpec("turbulent_prandtl",ParameterType.REAL,None,0.9),
    FieldSpec("k_min",ParameterType.REAL,"m^2/s^2",1.0e-12),
    FieldSpec("epsilon_min",ParameterType.REAL,"m^2/s^3",1.0e-12),
    FieldSpec("omega_min",ParameterType.REAL,"1/s",1.0e-12),
),("incompressible",))
PHYSICS_SPECS = (INCOMPRESSIBLE, ENERGY, TURBULENCE)
TURBULENCE_MODELS = ("LAMINAR","KEPSILON","SST","SMAGORINSKY","DES")

def default_parameters(spec: PhysicsSpec) -> tuple[Parameter, ...]:
    choices = TURBULENCE_MODELS if spec.key == "turbulence" else ()
    return tuple(Parameter(f.name,f.default,f.kind,f.unit,choices) for f in spec.fields)

def physics_spec(key: str) -> PhysicsSpec:
    for spec in PHYSICS_SPECS:
        if spec.key == key:
            return spec
    raise KeyError(key)
