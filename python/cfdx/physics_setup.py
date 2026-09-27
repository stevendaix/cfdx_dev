"""Structured setup schemas for supported CFDX physics controls.

This is an application/UI schema, not a second solver implementation. The
turbulence catalogue deliberately distinguishes solver-ready transport models
from algebraic/kernel-only models so the GUI never presents an unvalidated
kernel as a production solver option.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

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


class TurbulenceFamily(str, Enum):
    LAMINAR = "laminar"
    RANS = "rans"
    LES = "les"
    HYBRID = "hybrid"


class TurbulenceStatus(str, Enum):
    SOLVER_READY = "solver_ready"
    KERNEL_ONLY = "kernel_only"
    PLANNED = "planned"


@dataclass(frozen=True)
class TurbulenceModelSpec:
    key: str
    label: str
    family: TurbulenceFamily
    status: TurbulenceStatus
    equations: int | None
    wall_bounded: bool
    steady: bool
    transient: bool
    dimensions: tuple[int, ...]
    description: str
    caveats: tuple[str, ...] = ()


# This catalogue is intentionally richer than the current dropdown. It is the
# contract used by the GUI to explain why an option is shown, hidden, or marked
# experimental. Numerical evaluation remains in C++.
TURBULENCE_CATALOG: tuple[TurbulenceModelSpec, ...] = (
    TurbulenceModelSpec(
        "LAMINAR", "Laminar", TurbulenceFamily.LAMINAR,
        TurbulenceStatus.SOLVER_READY, 0, False, True, True, (2, 3),
        "No turbulence closure.",
    ),
    TurbulenceModelSpec(
        "KEPSILON", "k-epsilon (standard)", TurbulenceFamily.RANS,
        TurbulenceStatus.SOLVER_READY, 2, True, True, True, (2, 3),
        "Classical two-equation eddy-viscosity RANS closure.",
    ),
    TurbulenceModelSpec(
        "RNG_KEPSILON", "RNG k-epsilon", TurbulenceFamily.RANS,
        TurbulenceStatus.SOLVER_READY, 2, True, True, True, (2, 3),
        "RNG k-epsilon transport model with its dedicated production correction.",
    ),
    TurbulenceModelSpec(
        "REALIZABLE_KEPSILON", "Realizable k-epsilon", TurbulenceFamily.RANS,
        TurbulenceStatus.SOLVER_READY, 2, True, True, True, (2, 3),
        "Variable-Cmu realizable k-epsilon transport model.",
    ),
    TurbulenceModelSpec(
        "KOMEGA", "k-omega", TurbulenceFamily.RANS,
        TurbulenceStatus.SOLVER_READY, 2, True, True, True, (2, 3),
        "Two-equation k-omega transport model.",
    ),
    TurbulenceModelSpec(
        "SST", "k-omega SST", TurbulenceFamily.RANS,
        TurbulenceStatus.SOLVER_READY, 2, True, True, True, (2, 3),
        "Menter SST transport model with blending, cross-diffusion and production limiting.",
    ),
    TurbulenceModelSpec(
        "SPALART_ALLMARAS", "Spalart-Allmaras", TurbulenceFamily.RANS,
        TurbulenceStatus.SOLVER_READY, 1, True, True, True, (2, 3),
        "One-equation eddy-viscosity transport model.",
        ("SA-negative, rotation, compressibility and QCR variants are not yet exposed.",),
    ),
    TurbulenceModelSpec(
        "SMAGORINSKY", "Smagorinsky LES", TurbulenceFamily.LES,
        TurbulenceStatus.KERNEL_ONLY, 0, True, False, True, (3,),
        "Algebraic subgrid-scale eddy-viscosity kernel.",
        ("Not yet a complete transported LES solver.",),
    ),
    TurbulenceModelSpec(
        "WALE", "WALE LES", TurbulenceFamily.LES,
        TurbulenceStatus.KERNEL_ONLY, 0, True, False, True, (3,),
        "WALE algebraic subgrid-scale kernel.",
        ("Not yet a complete transported LES solver.",),
    ),
    TurbulenceModelSpec(
        "DYNAMIC_KEQN", "Dynamic one-equation LES", TurbulenceFamily.LES,
        TurbulenceStatus.PLANNED, 1, True, False, True, (3,),
        "Dynamic one-equation subgrid-scale model.",
        ("The enum exists, but a complete transport implementation is not currently exposed.",),
    ),
    TurbulenceModelSpec(
        "DES", "DES", TurbulenceFamily.HYBRID,
        TurbulenceStatus.KERNEL_ONLY, 0, True, False, True, (3,),
        "Algebraic DES length-scale/eddy-viscosity kernel.",
        ("Not yet a complete hybrid RANS-LES transport implementation.",),
    ),
    TurbulenceModelSpec(
        "DDES", "DDES", TurbulenceFamily.HYBRID,
        TurbulenceStatus.KERNEL_ONLY, 0, True, False, True, (3,),
        "DDES shielding and length-scale kernel.",
        ("Transport, shielding and wall-treatment integration are not complete.",),
    ),
    TurbulenceModelSpec(
        "IDDES", "IDDES", TurbulenceFamily.HYBRID,
        TurbulenceStatus.KERNEL_ONLY, 0, True, False, True, (3,),
        "IDDES shielding/blending length-scale kernel.",
        ("Transport, shielding and wall-treatment integration are not complete.",),
    ),
)

TURBULENCE_MODELS = tuple(
    model.key for model in TURBULENCE_CATALOG
    if model.status is TurbulenceStatus.SOLVER_READY
)


@dataclass(frozen=True)
class TurbulenceSelectionContext:
    """User-facing context used to produce transparent model suggestions.

    Suggestions are advisory only: the caller must explicitly select the
    resulting model. No numerical model is silently changed.
    """

    application: str = "general"
    dimensions: int = 3
    steady: bool = True
    wall_bounded: bool = True
    target_y_plus: float | None = None
    strong_separation: bool = False
    strong_rotation: bool = False
    transition_expected: bool = False
    high_fidelity: bool = False


@dataclass(frozen=True)
class TurbulenceRecommendation:
    model: str
    reasons: tuple[str, ...]
    warnings: tuple[str, ...] = ()


def turbulence_model(key: str) -> TurbulenceModelSpec:
    for model in TURBULENCE_CATALOG:
        if model.key == key:
            return model
    raise KeyError(key)


def recommend_turbulence_models(
    context: TurbulenceSelectionContext,
) -> tuple[TurbulenceRecommendation, ...]:
    """Return transparent, rule-based suggestions from the supported set.

    The rules intentionally prefer documented applicability over an opaque
    automatic selector. Future models can be added to the catalogue without
    changing the case-file format.
    """
    if context.dimensions not in (2, 3):
        raise ValueError("dimensions must be 2 or 3")
    if context.target_y_plus is not None and context.target_y_plus < 0:
        raise ValueError("target_y_plus must be non-negative")

    candidates = ["SST", "REALIZABLE_KEPSILON", "KOMEGA", "KEPSILON", "RNG_KEPSILON", "SPALART_ALLMARAS"]
    recommendations: list[TurbulenceRecommendation] = []

    for key in candidates:
        reasons: list[str] = []
        warnings: list[str] = []

        if key == "SST":
            reasons.append("general wall-bounded RANS choice with near-wall and adverse-pressure-gradient coverage")
            if context.strong_separation:
                reasons.append("SST is explicitly retained as a separation-oriented RANS candidate")
            if context.strong_rotation:
                warnings.append("a dedicated rotation/curvature correction is not yet exposed in CFDX")
        elif key == "SPALART_ALLMARAS":
            reasons.append("compact one-equation RANS model suited to attached external aerodynamic boundary layers")
            if context.application == "external":
                reasons.append("external-aerodynamics use is a natural target for SA")
            if context.strong_separation:
                warnings.append("separation behaviour should be qualified on the actual target geometry")
        elif key == "REALIZABLE_KEPSILON":
            reasons.append("two-equation RANS option for general industrial/internal-flow use")
        elif key == "RNG_KEPSILON":
            reasons.append("two-equation RANS variant with an RNG production correction")
        elif key == "KEPSILON":
            reasons.append("classical, economical two-equation RANS baseline")
        elif key == "KOMEGA":
            reasons.append("two-equation near-wall-oriented RANS alternative")

        if context.transition_expected:
            warnings.append("transition transport is not yet implemented; these are fully turbulent RANS models")
        if context.target_y_plus is not None and context.target_y_plus > 30:
            warnings.append("wall-function treatment is not yet a complete CFDX turbulence feature; do not infer wall-function readiness from y+ alone")
        if context.high_fidelity and context.dimensions == 3 and context.steady:
            warnings.append("high-fidelity scale-resolving models require transient 3-D solver support not yet complete in CFDX")
        if context.dimensions == 2 and key in {"SST", "SPALART_ALLMARAS", "REALIZABLE_KEPSILON", "KOMEGA", "KEPSILON", "RNG_KEPSILON"}:
            reasons.append("available in the current 2-D transport architecture")

        recommendations.append(TurbulenceRecommendation(key, tuple(reasons), tuple(warnings)))

    return tuple(recommendations)


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
    FieldSpec("model",ParameterType.CHOICE,None,"SST"),
    FieldSpec("density",ParameterType.REAL,"kg/m^3",1.0),
    FieldSpec("molecular_viscosity",ParameterType.REAL,"Pa s",1.0e-3),
    FieldSpec("turbulent_prandtl",ParameterType.REAL,None,0.9),
    FieldSpec("k_min",ParameterType.REAL,"m^2/s^2",1.0e-12),
    FieldSpec("epsilon_min",ParameterType.REAL,"m^2/s^3",1.0e-12),
    FieldSpec("omega_min",ParameterType.REAL,"1/s",1.0e-12),
),("incompressible",))
PHYSICS_SPECS = (INCOMPRESSIBLE, ENERGY, TURBULENCE)


def default_parameters(spec: PhysicsSpec) -> tuple[Parameter, ...]:
    choices = TURBULENCE_MODELS if spec.key == "turbulence" else ()
    return tuple(Parameter(f.name,f.default,f.kind,f.unit,choices) for f in spec.fields)


def physics_spec(key: str) -> PhysicsSpec:
    for spec in PHYSICS_SPECS:
        if spec.key == key:
            return spec
    raise KeyError(key)

def turbulence_model_from_case(values: dict[str, object]) -> dict[str, object]:
    # Normalize legacy flat turbulence settings to the structured contract.
    model = str(values.get("model", "SST")).upper()
    spec = turbulence_model(model)
    family = spec.family.value.upper()
    return {
        **values,
        "family": family,
        "model": model,
        "wall_treatment": values.get("wall_treatment", "resolved" if spec.wall_bounded else "none"),
        "transition": values.get("transition", {"model": "none"}),
        "corrections": dict(values.get("corrections", {})),
    }

def validate_turbulence_selection(values: dict[str, object]) -> dict[str, object]:
    normalized = turbulence_model_from_case(values)
    model = turbulence_model(str(normalized["model"]))
    if model.status is not TurbulenceStatus.SOLVER_READY:
        raise ValueError(f"{model.label} is {model.status.value}; it is not a production solver model")
    if normalized["wall_treatment"] not in {"resolved", "wall_function", "all_y_plus", "none"}:
        raise ValueError("unsupported turbulence wall treatment")
    if not isinstance(normalized["corrections"], dict):
        raise ValueError("turbulence corrections must be a mapping")
    if not isinstance(normalized["transition"], dict):
        raise ValueError("turbulence transition must be a mapping")
    return normalized