"""Minimal Fluent-like BC choices lowered to the native CFDX CaseSetup schema.

These are *physical* boundary choices, not a replacement for the FVM operators.
Unsupported Fluent BCs fail closed, rather than masquerading as velocity-inlet
or zero-gradient. This module is GUI-independent and can be used from a CLI.
"""

from __future__ import annotations

from math import isfinite
from typing import Sequence

from .io.schema import BCType, BCValueType, BoundarySpec

# Display label -> canonical physical choice. Only BCs already understood by
# the incompressible production adapter are offered as runnable presets.
SIMPLE_BOUNDARY_CHOICES = (
    ("velocity-inlet", "Velocity Inlet", "Vitesse et direction imposees"),
    ("pressure-outlet", "Pressure Outlet", "Pression statique de sortie"),
    ("no-slip-wall", "Wall / No Slip", "Paroi fixe"),
    ("moving-wall", "Wall / Moving", "Vitesse de paroi imposee"),
    ("symmetry", "Symmetry", "Flux normal nul (patch SYMMETRY requis)"),
    ("outflow", "Outflow", "Sortie en gradient nul, sans pression imposee"),
    ("empty", "Empty (2D)", "Face extrudee sans transport normal"),
)

UNSUPPORTED_CHOICE_REASONS = {
    "pressure-inlet": (
        "pression totale / relation de Bernoulli non integree au solveur incompressible"
    ),
    "mass-flow-inlet": (
        "conversion de debit en flux de face non integree et non qualifiee"
    ),
    "pressure-far-field": (
        "condition caracteristique Fluent pour l'ecoulement compressible, "
        "non equivalente au farfield incompressible CFDX"
    ),
    "slip-wall": (
        "operateur de non-penetration et contraintes tangentielles non raccorde"
    ),
    "periodic": (
        "appariement des faces periodiques non raccorde a ce solveur de production"
    ),
    "interface": (
        "couplage conservatif des faces d'interface non raccorde "
        "a ce solveur de production"
    ),
}


def make_simple_boundary(
    patch_name: str,
    kind: str,
    *,
    velocity: Sequence[float] | None = None,
    pressure_pa: float | None = None,
    wall_temperature_k: float | None = None,
) -> BoundarySpec:
    """Build a physical BC and validate the *exact* fields used by runtime.

    Returns a native `BoundarySpec`, serializable with `.model_dump(mode='json')`.
    For SST, CH/T and radiation use advanced native scalar BC contracts; do not
    silently populate their required fields with fictional default values.
    """
    if not patch_name or not patch_name.strip():
        raise ValueError("un nom de patch non vide est requis")
    if kind in UNSUPPORTED_CHOICE_REASONS:
        raise ValueError(f"{kind} non disponible : {UNSUPPORTED_CHOICE_REASONS[kind]}")
    available = {name for name, _, _ in SIMPLE_BOUNDARY_CHOICES}
    if kind not in available:
        raise ValueError(f"condition limite inconnue : {kind}")

    if velocity is None:
        vector = None
    else:
        vector = list(velocity)
        if len(vector) != 3 or any(
            isinstance(x, bool) or not isinstance(x, (int, float)) or not isfinite(x)
            for x in vector
        ):
            raise ValueError(
                "la vitesse doit contenir 3 composantes reelles finies (m/s)"
            )
        vector = [float(x) for x in vector]

    if pressure_pa is not None and (
        isinstance(pressure_pa, bool) or not isfinite(pressure_pa)
    ):
        raise ValueError("la pression de sortie doit etre finie (Pa)")
    if wall_temperature_k is not None and (
        isinstance(wall_temperature_k, bool)
        or not isinstance(wall_temperature_k, (int, float))
        or not isfinite(wall_temperature_k)
        or wall_temperature_k <= 0
    ):
        raise ValueError(
            "la temperature de paroi doit etre en kelvins et strictement positive"
        )

    needs_velocity = kind in ("velocity-inlet", "moving-wall")
    if needs_velocity != (vector is not None):
        raise ValueError(
            "une vitesse de 3 composantes est requise uniquement "
            "pour Velocity Inlet et Moving Wall"
        )
    if (pressure_pa is not None) != (kind == "pressure-outlet"):
        raise ValueError(
            "la pression relative est requise uniquement pour Pressure Outlet"
        )
    if wall_temperature_k is not None and kind != "no-slip-wall":
        raise ValueError(
            "la temperature de paroi ne s'applique qu'a Wall / No Slip dans ce preset"
        )

    spec = BoundarySpec(patch_name=patch_name.strip())
    if kind == "velocity-inlet":
        spec.type = BCType.INLET
        spec.value_type = BCValueType.FIXED
        spec.velocity_vector = vector
        spec.source_type_name = "velocity-inlet"
    elif kind == "pressure-outlet":
        spec.type = BCType.PRESSURE_OUTLET
        spec.value_type = BCValueType.OUTLET_PRESSURE
        spec.pressure = float(pressure_pa)
        spec.back_pressure = float(pressure_pa)
        spec.source_type_name = "pressure-outlet"
    elif kind == "no-slip-wall":
        spec.type = BCType.WALL
        spec.value_type = BCValueType.WALL_NO_SLIP
        if wall_temperature_k is not None:
            spec.value_type = BCValueType.WALL_THERMAL
            spec.wall_temperature = float(wall_temperature_k)
        spec.source_type_name = "wall"
    elif kind == "moving-wall":
        spec.type = BCType.WALL
        spec.value_type = BCValueType.FIXED
        spec.velocity_vector = vector
        spec.source_type_name = "wall"
    elif kind == "symmetry":
        spec.type = BCType.SYMMETRY
        spec.value_type = BCValueType.ZERO_GRADIENT
        spec.source_type_name = "symmetry"
    elif kind == "outflow":
        spec.type = BCType.OUTLET
        spec.value_type = BCValueType.ZERO_GRADIENT
        spec.source_type_name = "outflow"
    elif kind == "empty":
        spec.type = BCType.EMPTY
        spec.value_type = BCValueType.ZERO_GRADIENT
        spec.source_type_name = "empty"
    return spec
