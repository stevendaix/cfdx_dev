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
