"""Scientific operations overview for Workbench (headless / read-only).

The dashboard intentionally separates a configured case, a successful setup
check, the solver convergence certificate and independent physical validation.
It never promotes console messages or a successful process exit to evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .application.state import ApplicationState
    from .application.run_center import RunCenterState


@dataclass(frozen=True, slots=True)
class OverviewFinding:
    key: str
    title: str
    detail: str
    level: str  # information, action, warning, error
    workspace: str  # one of the Workbench's stable workspace keys

    def to_dict(self) -> dict[str, str]:
        return {
            "key": self.key, "title": self.title, "detail": self.detail,
            "level": self.level, "workspace": self.workspace,
        }


@dataclass(frozen=True, slots=True)
class OperationsOverview:
    project: str
    project_path: str | None
    unsaved_changes: bool
    execution_state: str
    iteration: int
    physical_time: float
    configured_sections: int
    section_count: int
    check_state: str
    check_scope: str
    solver_evidence: str
    solver_evidence_sha256: tuple[str, str] | None
    mesh_catalog_loaded: bool
    results_frames: int
    findings: tuple[OverviewFinding, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "format": "CFDX-WORKBENCH-OVERVIEW", "schema_version": 1,
            "project": self.project, "project_path": self.project_path,
            "unsaved_changes": self.unsaved_changes,
            "execution_state": self.execution_state,
            "iteration": self.iteration, "physical_time": self.physical_time,
            "configuration": {"sections_recorded": self.configured_sections,
                              "sections_considered": self.section_count,
                              "not_a_readiness_certificate": True},
            "setup_check": {"state": self.check_state, "scope": self.check_scope},
            "solver_evidence": {"state": self.solver_evidence,
                                "sha256": list(self.solver_evidence_sha256)
                                if self.solver_evidence_sha256 else None,
                                "independent_physical_qualification": False},
            "mesh_catalog_loaded": self.mesh_catalog_loaded,
            "results_frames": self.results_frames,
            "findings": [row.to_dict() for row in self.findings],
        }


def build_operations_overview(
    state: "ApplicationState", *,
    run: "RunCenterState | None" = None,
    mesh_catalog_loaded: bool = False,
) -> OperationsOverview:
    """Compute an honest presentation snapshot, without touching Qt or files.

    Workflow COMPLETE means *recorded configuration or declared availability*,
    NOT scientific qualification. Only a run-scoped ConvergenceEvidence may
    yield a `converged` solver-evidence state; that still cannot certify V&V.
    """
    findings: list[OverviewFinding] = []
    workflow = {step.id: step for step in state.workflow if step.parent is not None}
    tracked = ("physics", "materials", "boundaries", "numerics", "solver")
    configured = sum(
        1 for name in tracked
        if name in workflow and str(getattr(workflow[name].status, "value", workflow[name].status)) == "complete"
    )
    destinations = {
        "physics": "setup", "materials": "setup", "boundaries": "setup",
        "numerics": "setup", "solver": "setup", "mesh": "mesh",
        "geometry": "mesh", "execution": "run",
    }
    if state.project.path is None:
        findings.append(OverviewFinding(
            "unsaved", "Project not saved", "Save the setup before a native case-based run.",
            "action", "setup",
        ))
    if state.project.dirty:
        findings.append(OverviewFinding(
            "dirty", "Unsaved configuration changes",
            "The live case has changes not yet persisted to the project file.",
            "warning", "setup",
        ))
    if state.requires_restart:
        findings.append(OverviewFinding(
            "restart", "Restart required",
            "A change to the numerical setup requires restarting the solver.",
            "warning", "run",
        ))
    for name in tracked:
        step = workflow.get(name)
        if step is None:
            continue
        status = str(getattr(step.status, "value", step.status))
        if status in ("error", "warning", "not_configured", "conditional"):
            # Workflow messages describe known configuration / adapter limits.
            severity = {
                "error": "error", "warning": "warning",
                "not_configured": "action", "conditional": "information",
            }[status]
            findings.append(OverviewFinding(
                f"workflow.{name}", f"{step.label}: {status.replace('_', ' ')}",
                step.message or "Review this section; a recorded value is not a physics check.",
                severity, destinations[name],
            ))
    if not mesh_catalog_loaded:
        findings.append(OverviewFinding(
            "mesh.not_inspected", "Mesh catalog not loaded",
            "Open the case mesh to inspect patches and cell counts. A missing preview does not prove the mesh is invalid.",
            "information", "mesh",
        ))
    errors = sum(
        str(getattr(d, "severity", "")).lower() in ("error", "critical")
        for d in state.diagnostics
    )
    warnings = sum(str(getattr(d, "severity", "")).lower() == "warning" for d in state.diagnostics)
    if not state.validation_current:
        check_state = "stale" if state.validation_performed else "not_checked"
        findings.append(OverviewFinding(
            "checks", "Setup checks not current",
            "Run Check against the current configuration. An empty diagnostics list is not proof of validation.",
            "action", "setup",
        ))
    elif errors:
        check_state = "errors"
        findings.insert(0, OverviewFinding(
            "checks", f"{errors} setup error(s)",
            "Correct all reported setup errors before starting a production simulation.",
            "error", "setup",
        ))
    elif warnings:
        check_state = "warnings"
        findings.append(OverviewFinding(
            "checks", f"{warnings} setup warning(s)",
            "Review the setup check warnings; they have not been suppressed.",
            "warning", "setup",
        ))
    else:
        check_state = "passed"
    check_scope = "case+mesh" if state.validation_with_mesh else "case_only"

    evidence = None if run is None else run.convergence_evidence
    signatures = None
    if evidence is not None and getattr(run, "convergence_verdict", None) not in (
        None, evidence.overall_converged
    ):
        findings.append(OverviewFinding(
            "solver.contradiction", "Inconsistent solver verdict",
            "The displayed run verdict disagrees with the native evidence; no convergence is certified.",
            "error", "run",
        ))
        evidence = None
    if evidence is not None:
        # ExecutionController owns run-scoped freshness checks. This viewer
        # only accepts its already parsed and verified evidence object.
        signatures = (evidence.execution_sha256, evidence.convergence_sha256)
        verdict = "converged" if evidence.overall_converged else "not_converged"
        if not evidence.overall_converged:
            findings.append(OverviewFinding(
                "solver.not_converged", "Solver did not converge",
                evidence.reason[:500], "warning", "run",
            ))
    else:
        verdict = "missing_or_unverified"
        if run is not None and run.convergence_evidence_error:
            findings.append(OverviewFinding(
                "solver.evidence", "Unverified solver evidence",
                str(run.convergence_evidence_error)[:500], "warning", "run",
            ))
        elif state.execution.iteration > 0:
            findings.append(OverviewFinding(
                "solver.evidence", "No verified solver certificate",
                "Iteration progress or a clean process exit is not proof of convergence.",
                "information", "run",
            ))

    # Keep a stable, bounded list: warnings/errors first, action items next.
    rank = {"error": 0, "warning": 1, "action": 2, "information": 3}
    findings.sort(key=lambda item: (rank[item.level], item.key))
    frames = len(getattr(state.results, "frames", ()))
    return OperationsOverview(
        project=state.project.name,
        project_path=state.project.path,
        unsaved_changes=state.project.dirty,
        execution_state=str(getattr(state.simulation_state, "value", state.simulation_state)),
        iteration=state.execution.iteration,
        physical_time=float(state.execution.time),
        configured_sections=configured, section_count=len(tracked),
        check_state=check_state, check_scope=check_scope,
        solver_evidence=verdict, solver_evidence_sha256=signatures,
        mesh_catalog_loaded=bool(mesh_catalog_loaded),
        results_frames=frames, findings=tuple(findings[:32]),
    )
