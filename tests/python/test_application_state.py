from dataclasses import FrozenInstanceError

import pytest

from cfdx import CFDXSession, SimulationState
from cfdx.application import (
    ApplicationStateChanged,
    EventBus,
    ExecutionChanged,
    SelectionState,
    build_application_state,
)


def test_application_state_is_a_read_only_snapshot(tmp_path) -> None:
    session = CFDXSession()
    session.case.name = "channel"
    session.edit("cfl", 2.0)

    state = build_application_state(
        session,
        project_path=tmp_path / "channel.cfdx.h5",
        dirty=True,
        selection=SelectionState("patch:inlet", "patch", "inlet"),
        capabilities=["solver.steady", "solver.steady"],
    )

    assert state.project.name == "channel"
    assert state.project.path.endswith("channel.cfdx.h5")
    assert state.project.dirty is True
    assert state.numerics_revision == session.numerics_revision
    assert state.capabilities == ("solver.steady",)
    assert state.selection.label == "inlet"
    assert state.case_tree == session.case_tree()
    with pytest.raises(FrozenInstanceError):
        state.project = state.project


def test_application_state_reads_execution_controller_without_owning_it() -> None:
    session = CFDXSession()

    class Controller:
        latest_metrics = "metrics"
        error = None

        class Series:
            samples = [1, 2]

        monitor_series = Series()

    state = build_application_state(session, controller=Controller())
    assert state.execution.latest_metrics == "metrics"
    assert state.execution.monitor_sample_count == 2
    assert state.simulation_state is SimulationState.CREATED


def test_event_bus_delivers_specific_and_general_state_events() -> None:
    session = CFDXSession()
    state = build_application_state(session)
    bus = EventBus()
    received = []
    bus.subscribe(ExecutionChanged, received.append)
    bus.subscribe(ApplicationStateChanged, received.append)

    event = ExecutionChanged(state)
    bus.publish(event)

    assert received == [event, event]
    bus.unsubscribe(ExecutionChanged, received.append)
    bus.publish(ExecutionChanged(state))
    assert received == [event, event, event]
