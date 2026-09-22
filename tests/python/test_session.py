from __future__ import annotations

import pytest

from cfdx import CFDXSession, ChangeImpact, SimulationState, TuiRenderer


def test_session_lifecycle_and_checkpoint() -> None:
    session = CFDXSession()
    assert session.state is SimulationState.CREATED

    session.validate()
    assert session.state is SimulationState.READY

    session.run()
    session.advance(iterations=3, dt=0.1)
    assert session.iteration == 3
    assert session.time == pytest.approx(0.3)

    session.pause()
    assert session.state is SimulationState.PAUSED
    checkpoint = session.checkpoint()
    assert checkpoint.iteration == 3
    assert checkpoint.time == pytest.approx(0.3)


def test_hot_and_restart_edits_are_explicit() -> None:
    session = CFDXSession()
    session.validate()

    assert session.edit("output.frequency", 10) is ChangeImpact.HOT
    assert not session.requires_restart

    assert session.edit("pressure.solver", "gmres", ChangeImpact.REBUILD) is ChangeImpact.REBUILD
    assert session.requires_restart

    session.acknowledge_restart()
    assert not session.requires_restart
    assert session.numerics_revision == 1


def test_invalid_running_edit_is_rejected() -> None:
    session = CFDXSession()
    session.validate()
    session.run()

    with pytest.raises(RuntimeError, match="case edits"):
        session.edit("solver", "new")


def test_case_tree_is_stable() -> None:
    session = CFDXSession()
    ids = [node.id for node in session.case_tree()]
    assert ids == [
        "geometry", "mesh", "physics", "materials", "boundaries",
        "numerics", "solver", "run", "monitors", "results", "reports",
    ]


def test_tui_is_headless_and_reflects_state() -> None:
    session = CFDXSession()
    session.validate()
    session.run()
    session.pause()
    session.edit("solver", "new", ChangeImpact.RESTART)

    rendered = TuiRenderer.render(session)
    assert "CFDX | Case: untitled" in rendered
    assert "State: PAUSED" in rendered
    assert "[RUN]" in rendered
    assert "restart/rebuild required" in rendered
