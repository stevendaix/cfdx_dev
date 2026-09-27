"""Typed case setup/configuration schema for CFDX I/O converters.

This module defines the CFDX intermediate representation (IR) that all solver
adapters translate external formats into. The schema is versioned and uses
Pydantic for strong typing and validation.

Corresponds to the C++ headers:
  src/cfdx_io/case_schema.h
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field

CFDX_SCHEMA_VERSION: int = 1


class BCType(str, Enum):
    WALL = "wall"
    INLET = "inlet"
    OUTLET = "outlet"
    PRESSURE_OUTLET = "pressure_outlet"
    SYMMETRY = "symmetry"
    PERIODIC = "periodic"
    INTERFACE = "interface"
    EMPTY = "empty"
    INTERNAL = "internal"
    UNKNOWN = "unknown"


class BCValueType(str, Enum):
    FIXED = "fixed"
    ZERO_GRADIENT = "zero_gradient"
    MIXED = "mixed"
    OUTLET_PRESSURE = "outlet_pressure"
    WALL_NO_SLIP = "wall_no_slip"
    WALL_SLIP = "wall_slip"
    WALL_THERMAL = "wall_thermal"
    UNKNOWN = "unknown"


class Severity(str, Enum):
    SUPPORTED = "supported"
    APPROXIMATED = "approximated"
    UNSUPPORTED_NONBLOCK = "unsupported_nonblocking"
    UNSUPPORTED_BLOCK = "unsupported_blocking"
    UNAVAILABLE = "unavailable"


class SourceInfo(BaseModel):
    solver: str = ""
    version: str = ""
    case_path: str = ""
    case_name: str = ""
    format: str = ""


class MaterialSpec(BaseModel):
    name: str = "air"
    density: float = 1.225
    dynamic_viscosity: float = 1.789e-5
    thermal_conductivity: float = 0.024
    specific_heat: float = 1006.43
    molecular_weight: float = 28.97
    eos_model: str = "ideal_gas"
    turbulence_model: str = "laminar"
    extra_properties: dict[str, str] = Field(default_factory=dict)


class BoundarySpec(BaseModel):
    patch_name: str = ""
    type: BCType = BCType.UNKNOWN
    value_type: BCValueType = BCValueType.UNKNOWN
    velocity_magnitude: float = 0.0
    velocity_vector: list[float] = Field(default_factory=list)
    temperature: float = 0.0
    pressure: float = 0.0
    turbulence_intensity: float = 0.0
    turbulence_length_scale: float = 0.0
    roughness_height: float = 0.0
    wall_temperature: float = 0.0
    back_pressure: float = 0.0
    source_zone_name: str = ""
    source_type_name: str = ""
    raw_params: dict[str, str] = Field(default_factory=dict)

    @property
    def has_raw_unmapped(self) -> bool:
        return bool(self.raw_params)


class InitialCondition(BaseModel):
    velocity: float = 0.0
    pressure: float = 0.0
    temperature: float = 288.15
    velocity_vector: list[float] = Field(default_factory=list)
    scalar_fields: dict[str, float] = Field(default_factory=dict)


class NumericalScheme(BaseModel):
    momentum_scheme: str = "first_order"
    pressure_scheme: str = "standard"
    momentum_interpolation: str = "linear"
    transient_scheme: str = "steady"
    gradient_operator: str = "green_gauss_cell"
    under_relaxation_momentum: float = 0.7
    under_relaxation_pressure: float = 0.3
    coupled_solver: str = "SIMPLE"
    preconditioner: str = ""
    residual_target: str = "1e-5"
    max_iterations: int = 500
    raw_settings: dict[str, str] = Field(default_factory=dict)


class MeshMetadata(BaseModel):
    n_vertices: int = 0
    n_faces: int = 0
    n_cells: int = 0
    n_boundary_faces: int = 0
    n_internal_faces: int = 0
    n_patches: int = 0
    dimension: int = 3
    has_polyhedral: bool = False
    has_nonplanar_faces: bool = False
    mesh_type: str = "unstructured"
    cell_types: list[str] = Field(default_factory=list)


class CaseSetup(BaseModel):
    schema_version: int = Field(default=CFDX_SCHEMA_VERSION)
    source: SourceInfo = Field(default_factory=SourceInfo)
    mesh_info: MeshMetadata = Field(default_factory=MeshMetadata)

    # Physics
    physics_model: str = "incompressible_laminar"
    turbulence_model: str = "laminar"
    energy_model: str = "isothermal"
    multiphase_model: str = "none"
    radiation_model: str = "none"
    transient: bool = False
    time_step: float = 0.0
    end_time: float = 0.0
    max_time_steps: int = 0

    # Materials
    materials: list[MaterialSpec] = Field(default_factory=list)

    # Boundary conditions
    boundary_conditions: list[BoundarySpec] = Field(default_factory=list)

    # Initial conditions
    initial_condition: InitialCondition = Field(default_factory=InitialCondition)

    # Numerical schemes
    numerics: NumericalScheme = Field(default_factory=NumericalScheme)

    # Reference values
    ref_length: float = 1.0
    ref_density: float = 1.0
    ref_velocity: float = 1.0
    ref_temperature: float = 293.15
    ref_pressure: float = 101325.0

    # Other
    gravity_vector: str = ""
    units: str = "SI"
    solver_mode: str = "steady"
    source_metadata: dict[str, str] = Field(default_factory=dict)

    def find_boundary(self, name: str) -> Optional[BoundarySpec]:
        for bc in self.boundary_conditions:
            if bc.patch_name == name:
                return bc
        return None

    def find_material(self, name: str) -> Optional[MaterialSpec]:
        for m in self.materials:
            if m.name == name:
                return m
        return None
