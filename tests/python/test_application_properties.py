from cfdx import CFDXSession
from cfdx.application import Application, WorkflowStatus
from cfdx.setup_model import SetupDiagnostic


def test_selection_resolves_properties_and_edit_updates_authoritative_case() -> None:
    session = CFDXSession()
    session.case.numerics["cfl"] = 1.0
    application = Application(session)

    application.select("numerics", "numerics", "Numerics")
    properties = application.properties()
    assert [prop.key for prop in properties] == ["numerics.cfl"]

    state = application.set_property("numerics.cfl", 2.5)
    assert session.case.numerics["cfl"] == 2.5
    assert state.project.dirty is True
    assert state.selection.stable_id == "numerics"


def test_validation_diagnostics_map_paths_to_workflow_status() -> None:
    session = CFDXSession()
    application = Application(session)
    diagnostic = SetupDiagnostic("error", "CFL", "CFL must be positive", "numerics.cfl")

    state = application.set_diagnostics([diagnostic])
    numerics = next(step for step in state.workflow if step.id == "numerics")

    assert numerics.status is WorkflowStatus.ERROR
    assert numerics.message == "CFL must be positive"


def test_solver_property_marks_restart_required() -> None:
    session = CFDXSession()
    application = Application(session)

    application.select("solver", "solver", "Solver")
    state = application.set_property("execution.solver", "/tmp/solver")

    assert session.case.execution.solver == "/tmp/solver"
    assert state.requires_restart is True
