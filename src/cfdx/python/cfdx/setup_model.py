"""Headless setup contracts consumed by CLI, TUI and GUI.

The objects in this module describe user-editable values without introducing a
second numerical model. Numerical ownership stays with the C++ mesh/field and
physics implementations.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Any, Callable

class ChangeImpact(str, Enum):
    HOT = "Hot"
    REQUIRES_RESTART = "RequiresRestart"
    REQUIRES_REBUILD = "RequiresRebuild"

class ParameterType(str, Enum):
    BOOL = "bool"
    INTEGER = "integer"
    REAL = "real"
    STRING = "string"
    CHOICE = "choice"

Validator = Callable[[Any], str | None]

@dataclass(frozen=True)
class Parameter:
    """A typed user-facing parameter contract."""
    name: str
    value: Any
    kind: ParameterType
    unit: str | None = None
    choices: tuple[Any, ...] = ()
    validator: Validator | None = None
    impact: ChangeImpact = ChangeImpact.HOT

    def validate(self) -> None:
        if self.kind is ParameterType.BOOL and not isinstance(self.value, bool):
            raise ValueError(f"{self.name}: expected bool")
        if self.kind is ParameterType.INTEGER and (not isinstance(self.value, int) or isinstance(self.value, bool)):
            raise ValueError(f"{self.name}: expected integer")
        if self.kind is ParameterType.REAL and (isinstance(self.value, bool) or not isinstance(self.value, (int, float)) or not math.isfinite(float(self.value))):
            raise ValueError(f"{self.name}: expected finite real")
        if self.kind in {ParameterType.STRING, ParameterType.CHOICE} and not isinstance(self.value, str):
            raise ValueError(f"{self.name}: expected string")
        if self.choices and self.value not in self.choices:
            raise ValueError(f"{self.name}: unsupported value {self.value!r}")
        if self.validator is not None:
            message = self.validator(self.value)
            if message:
                raise ValueError(f"{self.name}: {message}")

@dataclass(frozen=True)
class MeshSelection:
    """Stable application identity for a mesh object."""
    kind: str
    index: int
    name: str | None = None
    stable_id: str | None = None

    def __post_init__(self) -> None:
        if self.kind not in {"body", "face", "edge", "region", "patch", "cell"}:
            raise ValueError(f"unsupported mesh selection kind: {self.kind}")
        if self.index < 0:
            raise ValueError("mesh selection index must be non-negative")
        if self.stable_id is not None and not self.stable_id.strip():
            raise ValueError("mesh selection stable_id must not be empty")

@dataclass(frozen=True)
class SetupDiagnostic:
    severity: str
    code: str
    message: str
    path: str = ""

    def __post_init__(self) -> None:
        if self.severity not in {"error", "warning"}:
            raise ValueError("severity must be error or warning")

def typed_parameter(name: str, value: Any, *, unit: str | None = None, choices: tuple[Any, ...] = (), impact: ChangeImpact = ChangeImpact.HOT) -> Parameter:
    """Infer the supported user-facing parameter type."""
    if isinstance(value, bool):
        kind = ParameterType.BOOL
    elif isinstance(value, int):
        kind = ParameterType.INTEGER
    elif isinstance(value, float):
        kind = ParameterType.REAL
    elif isinstance(value, str):
        kind = ParameterType.CHOICE if choices else ParameterType.STRING
    else:
        raise TypeError(f"unsupported parameter value type: {type(value).__name__}")
    parameter = Parameter(name, value, kind, unit, choices, None, impact)
    parameter.validate()
    return parameter
