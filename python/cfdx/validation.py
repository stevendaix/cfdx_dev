"""Headless case validation shared by CLI, TUI and GUI."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Protocol

from .case import Case
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

def validate_case(case: Case, mesh: MeshInfo | None = None) -> ValidationReport:
    diagnostics: list[SetupDiagnostic] = []
    if not case.name.strip():
        diagnostics.append(SetupDiagnostic("error", "CASE_NAME", "case name is empty", "name"))
    ranks = case.execution.mpi_ranks
    if not isinstance(ranks, int) or ranks < 1:
        diagnostics.append(SetupDiagnostic("error", "MPI_RANKS", "mpi_ranks must be a positive integer", "execution.mpi_ranks"))
    cfl = case.numerics.get("cfl")
    if cfl is not None and (isinstance(cfl, bool) or not isinstance(cfl, (int, float)) or not isfinite(float(cfl)) or cfl <= 0):
        diagnostics.append(SetupDiagnostic("error", "CFL", "CFL must be a finite positive number", "numerics.cfl"))
    if not case.execution.solver:
        diagnostics.append(SetupDiagnostic("warning", "SOLVER_UNSET", "no solver executable is configured", "execution.solver"))
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
    if mesh is not None:
        if mesh.n_cells <= 0:
            diagnostics.append(SetupDiagnostic("error", "EMPTY_MESH", "mesh contains no cells", "mesh"))
        mesh_patches = set(mesh.patches)
        for name in case.boundaries:
            if name not in mesh_patches:
                diagnostics.append(SetupDiagnostic("error", "ORPHAN_BOUNDARY", f"boundary {name!r} is not present in the loaded mesh", f"boundaries.{name}"))
    return ValidationReport(tuple(diagnostics))
