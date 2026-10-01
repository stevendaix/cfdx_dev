import pytest

from cfdx import CFDXSession, ChangeImpact, SimulationState
from cfdx.application import Application, ApplicationStateChanged


class FakeController:
    def __init__(self, session):
        self.session = session
        self.calls = []

    def start(self):
        self.calls.append("start")
        self.session.run()

    def pause(self):
        self.calls.append("pause")
        self.session.pause()

    def resume(self):
        self.calls.append("resume")
        self.session.run()

    def stop(self):
        self.calls.append("stop")
        self.session.stop()


def test_application_routes_configuration_through_a_command() -> None:
    session = CFDXSession()
    application = Application(session)

    state = application.set_numerical_option("cfl", 3.5, ChangeImpact.RESTART)

    assert session.case.numerics["cfl"] == 3.5
    assert state.dirty is True
    assert state.requires_restart is True


def test_application_routes_execution_to_controller_and_publishes_state() -> None:
    session = CFDXSession()
    controller = FakeController(session)
    application = Application(session, controller=controller)
    states = []
    application.events.subscribe(ApplicationStateChanged, states.append)

    application.run()
    application.pause()
    application.resume()
    application.stop()

    assert controller.calls == ["start", "pause", "resume", "stop"]
    assert session.state is SimulationState.STOPPED
    assert application.state.simulation_state is SimulationState.STOPPED
    assert len(states) == 4
    assert all(event.state.simulation_state is session.state for event in states)


def test_execution_command_requires_controller() -> None:
    with pytest.raises(RuntimeError, match="execution controller"):
        Application(CFDXSession()).run()
