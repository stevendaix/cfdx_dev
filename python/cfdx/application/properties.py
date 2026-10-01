"""Contextual property contracts for the Workbench Properties panel."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..case import Case
from ..setup_model import ChangeImpact, Parameter, ParameterType


@dataclass(frozen=True)
class PropertyState:
    key: str
    label: str
    value: Any
    kind: ParameterType
    unit: str | None = None
    editable: bool = True
    impact: ChangeImpact = ChangeImpact.HOT
    message: str | None = None


def properties_for_selection(case: Case, stable_id: str | None) -> tuple[PropertyState, ...]:
    """Resolve a stable navigation ID to editable application properties."""
    if stable_id == "numerics":
        value = case.numerics.get("cfl", 1.0)
        return (PropertyState("numerics.cfl", "CFL", value, ParameterType.REAL, impact=ChangeImpact.HOT),)
    if stable_id == "solver":
        return (PropertyState("execution.solver", "Solver executable", case.execution.solver or "", ParameterType.STRING, impact=ChangeImpact.REQUIRES_RESTART),)
    if stable_id == "physics":
        return tuple(
            PropertyState(f"physics.{name}.enabled", name.title(), bool(value.get("enabled", True)), ParameterType.BOOL)
            for name, value in case.physics.items() if isinstance(value, dict)
        )
    if stable_id == "materials":
        return tuple(
            PropertyState(f"materials.{name}", name.title(), value, ParameterType.STRING, editable=False)
            for name, value in case.materials.items()
        )
    if stable_id == "boundaries":
        return tuple(
            PropertyState(f"boundaries.{name}.type", name, value.get("type", ""), ParameterType.STRING, editable=False)
            for name, value in case.boundaries.items() if isinstance(value, dict)
        )
    return ()


def set_property(case: Case, key: str, value: Any) -> ChangeImpact:
    """Apply an explicitly routed property to the authoritative Case."""
    parts = key.split(".")
    if key == "numerics.cfl":
        parameter = Parameter("cfl", value, ParameterType.REAL)
        parameter.validate()
        case.numerics["cfl"] = value
        return ChangeImpact.HOT
    if key == "execution.solver":
        if not isinstance(value, str):
            raise ValueError("solver executable must be a string")
        case.execution.solver = value or None
        return ChangeImpact.REQUIRES_RESTART
    if len(parts) == 3 and parts[0] == "physics" and parts[2] == "enabled":
        if not isinstance(value, bool):
            raise ValueError("physics enabled must be boolean")
        model = case.physics.setdefault(parts[1], {})
        if not isinstance(model, dict):
            raise ValueError(f"physics.{parts[1]} must be a mapping")
        model["enabled"] = value
        return ChangeImpact.REQUIRES_RESTART
    raise KeyError(f"property is not editable: {key}")
