"""Initialization capabilities exposed by the application layer."""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from .setup_model import ChangeImpact

class InitializationMode(str, Enum):
    UNIFORM = "uniform"
    FIELD = "field"

@dataclass(frozen=True)
class InitializationSpec:
    mode: InitializationMode
    field: str | None = None
    value: float | None = None

    @property
    def change_impact(self) -> ChangeImpact:
        return ChangeImpact.REQUIRES_RESTART

    def validate(self) -> None:
        if self.mode is InitializationMode.UNIFORM:
            if self.value is None:
                raise ValueError("uniform initialization requires a value")
            if self.field is None or not self.field.strip():
                raise ValueError("uniform initialization requires a field")
        elif self.mode is InitializationMode.FIELD:
            if self.field is None or not self.field.strip():
                raise ValueError("field initialization requires a source field")

def supported_initialization_modes() -> tuple[InitializationMode,...]:
    # Standard/hybrid initialization is intentionally absent until the solver
    # exposes a corresponding numerical contract.
    return (InitializationMode.UNIFORM, InitializationMode.FIELD)
