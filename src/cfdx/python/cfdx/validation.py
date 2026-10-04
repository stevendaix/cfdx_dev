"""Headless case validation shared by CLI, TUI and GUI."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Protocol

from .boundary_setup import validate_boundary_definition
from .case import Case
from .initialization import InitializationMode, InitializationSpec
from .materials import MaterialSpec
from .physics_setup import PHYSICS_SPECS, physics_spec
from .setup_model import SetupDiagnostic


class MeshInfo(Protocol):
    @property
    def n_cells(self) -> int: ...

    @property
    def patches(self) -> tuple[str, ...]: ...


@dataclass(frozen=True)
class ValidationReport:
    diagnostics: tuple[SetupDiagnostic, ...]

    @property
    def errors(self) -> tuple[SetupDiagnostic, ...]:
        return tuple(d for d in self.diagnostics if d.severity == "error")

    @property
    def warnings(self) -> tuple[SetupDiagnostic, ...]:
        return tuple(d for d in self.diagnostics if d.severity == "warning")

    @property
    def ok(self) -> bool:
        return not self.errors


def _enabled_models(case: Case) -> set[str]:
    return {
        key for key, value in case.physics.items()
        if isinstance(value, dict) and bool(value.get("enabled", True))
    }


def validate_case(case: Case, mesh: MeshInfo | None = None) -> ValidationReport:
    diagnostics: list[SetupDiagnostic] = []

    if not case.name.strip():
        diagnostics.append(SetupDiagnostic("error", "CASE_NAME", "case name is empty", "name"))

    ranks = case.execution.mpi_ranks
    if not isinstance(ranks, int) or ranks < 1:
        diagnostics.append(SetupDiagnostic("error", "MPI_RANKS", "mpi_ranks must be a positive integer", "execution.mpi_ranks"))

    cfl = case.numerics.get("cfl")
    if cfl is not None and (
        isinstance(cfl, bool)
        or not isinstance(cfl, (int, float))
        or not isfinite(float(cfl))
        or cfl <= 0
    ):
        diagnostics.append(SetupDiagnostic("error", "CFL", "CFL must be a finite positive number", "numerics.cfl"))

    if not case.execution.solver:
        diagnostics.append(SetupDiagnostic("warning", "SOLVER_UNSET", "no solver executable is configured", "execution.solver"))

    if mesh is not None and mesh.n_cells <= 0:
        diagnostics.append(SetupDiagnostic("error", "EMPTY_MESH", "mesh contains no cells", "mesh"))

    enabled = _enabled_models(case)
    if not enabled:
        diagnostics.append(SetupDiagnostic("error", "PHYSICS_UNSET", "no enabled physics model is configured", "physics"))

    known_models = {spec.key: spec for spec in PHYSICS_SPECS}
    for model, value in case.physics.items():
        if model == "initialization":
            continue
        if model not in known_models:
            diagnostics.append(SetupDiagnostic("error", "PHYSICS_UNSUPPORTED", f"unsupported physics model {model!r}", f"physics.{model}"))
            continue
        if not isinstance(value, dict):
            diagnostics.append(SetupDiagnostic("error", "PHYSICS_SCHEMA", "physics definition must be a mapping", f"physics.{model}"))
            continue
        if not bool(value.get("enabled", True)):
            continue
        spec = known_models[model]
        missing = [dependency for dependency in spec.requires if dependency not in enabled]
        for dependency in missing:
            diagnostics.append(
                SetupDiagnostic(
                    "error",
                    "PHYSICS_DEPENDENCY",
                    f"{model} requires enabled physics model {dependency}",
                    f"physics.{model}",
                )
            )
        supported_fields = {field.name: field for field in spec.fields}
        for name, parameter in value.items():
            if name == "enabled":
                continue
            field = supported_fields.get(name)
            if field is None:
                diagnostics.append(
                    SetupDiagnostic("error", "PHYSICS_PARAMETER", f"unsupported parameter {name!r} for {model}", f"physics.{model}.{name}")
                )
                continue
            try:
                from .setup_model import Parameter
                Parameter(name, parameter, field.kind, field.unit).validate()
            except (TypeError, ValueError) as exc:
                diagnostics.append(SetupDiagnostic("error", "PHYSICS_VALUE", str(exc), f"physics.{model}.{name}"))

    for name, raw_material in case.materials.items():
        if not isinstance(raw_material, dict):
            diagnostics.append(SetupDiagnostic("error", "MATERIAL_SCHEMA", "material definition must be a mapping", f"materials.{name}"))
            continue
        try:
            MaterialSpec(
                name,
                float(raw_material.get("density", 0.0)),
                float(raw_material.get("dynamic_viscosity", 0.0)),
                float(raw_material.get("cp", 0.0)),
                float(raw_material.get("conductivity", 0.0)),
            ).validate()
        except (TypeError, ValueError) as exc:
            diagnostics.append(SetupDiagnostic("error", "MATERIAL_VALUE", str(exc), f"materials.{name}"))

    if enabled & {"incompressible", "energy", "turbulence"} and not case.materials:
        diagnostics.append(SetupDiagnostic("error", "MATERIAL_UNSET", "enabled physics requires at least one material", "materials"))

    initialization = case.physics.get("initialization")
    if initialization is not None:
        if not isinstance(initialization, dict):
            diagnostics.append(SetupDiagnostic("error", "INITIALIZATION_SCHEMA", "initialization definition must be a mapping", "physics.initialization"))
        else:
            try:
                InitializationSpec(
                    InitializationMode(initialization["mode"]),
                    initialization.get("field"),
                    initialization.get("value"),
                ).validate()
            except (KeyError, TypeError, ValueError) as exc:
                diagnostics.append(SetupDiagnostic("error", "INITIALIZATION", str(exc), "physics.initialization"))

    allowed = {"inlet", "outlet", "wall", "symmetry", "periodic", "interface", "empty"}
    for name, boundary in case.boundaries.items():
        if not name.strip():
            diagnostics.append(SetupDiagnostic("error", "BOUNDARY_NAME", "boundary name is empty", "boundaries"))
        if not isinstance(boundary, dict):
            diagnostics.append(SetupDiagnostic("error", "BOUNDARY_TYPE", "boundary definition must be a mapping", f"boundaries.{name}"))
            continue
        boundary_type = boundary.get("type")
        if boundary_type not in allowed:
            diagnostics.append(SetupDiagnostic("error", "BOUNDARY_KIND", f"unsupported boundary type {boundary_type!r}", f"boundaries.{name}.type"))

    for name, boundary in case.boundaries.items():
        try:
            validate_boundary_definition(
                boundary,
                known_fields=tuple(boundary.get("fields", ())) if isinstance(boundary, dict) else (),
            )
        except (KeyError, ValueError) as exc:
            diagnostics.append(SetupDiagnostic("error", "BOUNDARY_SCHEMA", str(exc), f"boundaries.{name}"))

    if mesh is not None:
        mesh_patches = {patch.name if hasattr(patch, "name") else patch for patch in mesh.patches}
        for name in case.boundaries:
            if name not in mesh_patches:
                diagnostics.append(
                    SetupDiagnostic(
                        "error",
                        "ORPHAN_BOUNDARY",
                        f"boundary {name!r} is not present in the loaded mesh",
                        f"boundaries.{name}",
                    )
                )

    return ValidationReport(tuple(diagnostics))
