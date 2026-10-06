"""Schema-driven setup contracts shared by GUI, TUI and CLI.

This module contains metadata only. It does not introduce a second CFD model:
the authoritative values remain on Case and existing physics/boundary schemas
remain the source of supported capabilities.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from .boundary_setup import fields_for_boundary_type
from .case import Case
from .physics_setup import physics_spec
from .setup_model import ChangeImpact, ParameterType

ValueReader = Callable[[Case], Any]

@dataclass(frozen=True)
class SetupField:
    """Stable metadata for one user-facing setup property."""
    key: str
    label: str
    kind: ParameterType
    reader: ValueReader
    unit: str | None = None
    choices: tuple[Any, ...] = ()
    editable: bool = True
    impact: ChangeImpact = ChangeImpact.HOT

    def value(self, case: Case) -> Any:
        return self.reader(case)

@dataclass(frozen=True)
class SetupDomain:
    """A named group of setup fields exposed by the application contract."""
    key: str
    label: str
    fields: tuple[SetupField, ...]

@dataclass(frozen=True)
class SetupSchema:
    """Immutable schema consumed by generic setup editors."""
    domains: tuple[SetupDomain, ...]

    def domain(self, key: str) -> SetupDomain:
        for domain in self.domains:
            if domain.key == key:
                return domain
        raise KeyError(key)

    def field(self, key: str) -> SetupField:
        for domain in self.domains:
            for field in domain.fields:
                if field.key == key:
                    return field
        raise KeyError(key)

def _physics_fields(case: Case) -> tuple[SetupField, ...]:
    fields: list[SetupField] = []
    for spec_key in ("incompressible", "energy", "turbulence"):
        spec = physics_spec(spec_key)
        model = case.physics.get(spec_key)
        if not isinstance(model, dict) or not model.get("enabled", False):
            continue
        for field in spec.fields:
            key = f"physics.{spec_key}.{field.name}"
            fields.append(
                SetupField(
                    key,
                    f"{spec.label}: {field.name.replace('_', ' ').title()}",
                    field.kind,
                    lambda c, sk=spec_key, fn=field.name: c.physics[sk].get(fn, field.default),
                    field.unit,
                    (),
                    True,
                    ChangeImpact.REQUIRES_RESTART,
                )
            )
    return tuple(fields)

def _boundary_fields(case: Case) -> tuple[SetupField, ...]:
    fields: list[SetupField] = []
    for name, boundary in sorted(case.boundaries.items()):
        if not isinstance(boundary, dict):
            continue
        boundary_type = boundary.get("type")
        if not isinstance(boundary_type, str):
            continue
        for spec in fields_for_boundary_type(boundary_type):
            if spec.name == "type":
                continue
            key = f"boundaries.{name}.{spec.name}"
            fields.append(
                SetupField(
                    key,
                    f"{name}: {spec.name.replace('_', ' ').title()}",
                    spec.kind,
                    lambda c, bn=name, fn=spec.name: c.boundaries[bn].get(fn),
                    spec.unit,
                    (),
                    True,
                    spec.impact,
                )
            )
    return tuple(fields)

def build_setup_schema(case: Case) -> SetupSchema:
    """Build the current schema without copying or mutating case state."""
    numerics = SetupDomain(
        "numerics",
        "Numerics",
        (SetupField("numerics.cfl", "CFL", ParameterType.REAL,
                    lambda c: c.numerics.get("cfl", 1.0), impact=ChangeImpact.HOT),),
    )
    solver = SetupDomain(
        "solver",
        "Solver",
        (SetupField("execution.solver", "Solver executable", ParameterType.STRING,
                    lambda c: c.execution.solver or "", impact=ChangeImpact.REQUIRES_RESTART),),
    )
    return SetupSchema((
        SetupDomain("physics", "Physics", _physics_fields(case)),
        SetupDomain("boundaries", "Boundary conditions", _boundary_fields(case)),
        numerics,
        solver,
    ))

def setup_fields_for_selection(case: Case, stable_id: str | None) -> tuple[SetupField, ...]:
    """Resolve a navigation domain to fields of the shared setup schema."""
    if stable_id is None:
        return ()
    return build_setup_schema(case).domain(stable_id).fields
