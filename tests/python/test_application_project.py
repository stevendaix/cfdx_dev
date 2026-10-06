from pathlib import Path

from cfdx import CFDXSession
from cfdx.application.application import Application


def test_application_project_lifecycle_uses_shared_case_boundary(tmp_path: Path) -> None:
    path = tmp_path / "channel.cfdx.h5"
    session = CFDXSession()
    session.case.name = "channel"
    session.case.enable("incompressible")
    session.case.set_numerics(cfl=0.8)

    app = Application(session)
    saved = app.save_project(path)

    assert saved.project.path == str(path)
    assert saved.project.dirty is False

    opened = Application().open_project(path)
    assert opened.project.path == str(path)
    assert opened.project.name == "channel"
    assert opened.project.dirty is False
    assert opened.case_revision == session.case_revision
    assert opened.simulation_state.name == "CREATED"


def test_application_requires_project_path_for_save() -> None:
    app = Application(CFDXSession())

    try:
        app.save_project()
    except ValueError as exc:
        assert "project path" in str(exc)
    else:
        raise AssertionError("save_project() must require a path for a new application")
