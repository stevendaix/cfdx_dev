from pathlib import Path
from unittest.mock import patch

from cfdx import CFDXSession
from cfdx.application import Application, ApplicationStateChanged


class MinimalController:
    def __init__(self, session):
        self.session = session
        self.latest_metrics = None
        self.monitor_series = type("Series", (), {"samples": []})()
        self.error = None


def test_application_can_attach_controller_after_project_save() -> None:
    application = Application(CFDXSession())
    events = []
    application.events.subscribe(ApplicationStateChanged, events.append)

    state = application.attach_controller(MinimalController(application.session))

    assert application.run_center is not None
    assert state.execution.monitor_sample_count == 0
    assert events[-1].state.execution.monitor_sample_count == 0


def test_application_replaces_session_and_clears_project_context() -> None:
    application = Application(CFDXSession(), project_path="case.cfdx.h5")
    application.dirty = True

    loaded = CFDXSession()
    loaded.case.name = "loaded"
    state = application.replace_session(loaded, project_path="loaded.cfdx.h5")

    assert state.project.name == "loaded"
    assert state.project.path.endswith("loaded.cfdx.h5")
    assert state.project.dirty is False
    assert application.controller is None


def test_application_owns_controller_construction_and_restart_selection(
    tmp_path,
) -> None:
    application = Application(CFDXSession(), project_path=tmp_path / "case.cfdx.h5")
    application.session.case.execution.solver = "/usr/bin/cfdx-solver"
    application.session.case.execution.restart_option = "--restart"
    dat_path = tmp_path / "checkpoint.dat.h5"

    # The DAT parser is intentionally not bypassed by the application facade.
    # A valid artifact is not required to test the ownership boundary.
    class Restart:
        iteration = 12
        time = 0.25
        fields = ("p",)

    with patch("cfdx.application.application.read_dat_restart", return_value=Restart()):
        loaded = application.load_dat(dat_path, use_for_restart=True)

    assert loaded.iteration == 12
    assert application.restart_dat == dat_path

    controller = application.ensure_controller()
    assert controller.session is application.session
    assert controller.runner.command == (
        "/usr/bin/cfdx-solver",
        str(tmp_path / "case.cfdx.h5"),
        "--restart",
        str(dat_path),
    )


def test_application_new_project_clears_restart_and_controller() -> None:
    application = Application(CFDXSession(), project_path="case.cfdx.h5")
    application._restart_dat = Path("checkpoint.dat.h5")
    application.new_project()

    assert application.restart_dat is None
    assert application.controller is None


def test_application_open_project_with_dat_preserves_restart_identity(tmp_path):
    application = Application(CFDXSession())
    project_path = tmp_path / "case.cfdx.h5"
    dat_path = tmp_path / "case.dat.h5"

    loaded = CFDXSession()
    loaded.case.name = "loaded"

    class ProjectStub:
        path = project_path

        def load_with_dat(self, requested):
            assert requested == dat_path
            return loaded, dat_path

        def load(self):
            return loaded

    with patch("cfdx.application.application.Project", return_value=ProjectStub()):
        state = application.open_project(project_path, with_dat=True, dat_path=dat_path)

    assert state.project.path.endswith("case.cfdx.h5")
    assert application.restart_dat == dat_path
