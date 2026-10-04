"""Material property schema for professional CFDX setup."""
from __future__ import annotations
from dataclasses import dataclass
from .setup_model import Parameter, ParameterType

@dataclass(frozen=True)
class MaterialSpec:
    name: str
    density: float
    dynamic_viscosity: float
    cp: float = 1000.0
    conductivity: float = 0.0

    def parameters(self) -> tuple[Parameter,...]:
        return (
            Parameter("density",self.density,ParameterType.REAL,"kg/m^3"),
            Parameter("dynamic_viscosity",self.dynamic_viscosity,ParameterType.REAL,"Pa s"),
            Parameter("cp",self.cp,ParameterType.REAL,"J/(kg K)"),
            Parameter("conductivity",self.conductivity,ParameterType.REAL,"W/(m K)"),
        )

    def validate(self) -> None:
        if not self.name.strip(): raise ValueError("material name must not be empty")
        if self.density <= 0.0: raise ValueError("density must be positive")
        if self.dynamic_viscosity < 0.0: raise ValueError("dynamic viscosity must be non-negative")
        if self.cp <= 0.0: raise ValueError("specific heat must be positive")
        if self.conductivity < 0.0: raise ValueError("conductivity must be non-negative")
        for parameter in self.parameters(): parameter.validate()
