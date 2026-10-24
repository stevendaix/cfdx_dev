"""Behavior tests for the bounded application edit history.

Covers gap G4 of docs/gui/PHASE0_AUDIT.md: reversible edits recorded by the
application layer, bounded depth, and an undo/redo surface shared by every
frontend (GUI, CLI, TUI).
"""

import pytest
from cfdx.application.events import ApplicationStateChanged
from cfdx.application.history import EditHistory, HistoryEntry
from cfdx.application.properties import get_property

from cfdx import CFDXSession
from cfdx.application import Application


def _entry(label: str, log: list[str]) -> HistoryEntry:
    """Build one reversible entry appending its label to the shared log."""

    return HistoryEntry(
        description=label,
        undo=lambda: log.append(f"undo {label}"),
        redo=lambda: log.append(f"redo {label}"),
    )


def test_undo_reverts_entries_in_reverse_order() -> None:
    log: list[str] = []
    history = EditHistory()
    for label in ("first", "second", "third"):
        history.push(_entry(label, log))
    assert history.undo() == "third"
    assert history.undo() == "second"
    assert history.undo() == "first"
    assert log == ["undo third", "undo second", "undo first"]
    with pytest.raises(IndexError, match="nothing to undo"):
        history.undo()


def test_redo_replays_the_last_undone_entry() -> None:
    log: list[str] = []
    history = EditHistory()
    history.push(_entry("first", log))
    history.push(_entry("second", log))
    assert history.undo() == "second"
    assert history.can_redo
    assert history.redo() == "second"
    assert log == ["undo second", "redo second"]
    with pytest.raises(IndexError, match="nothing to redo"):
        history.redo()


def test_new_edit_clears_the_redo_branch() -> None:
    log: list[str] = []
    history = EditHistory()
    history.push(_entry("first", log))
    history.push(_entry("second", log))
    assert history.undo() == "second"
    history.push(_entry("third", log))
    assert not history.can_redo
    assert history.undo() == "third"
    assert history.undo() == "first"
    assert log == ["undo second", "undo third", "undo first"]


def test_undo_depth_stays_bounded() -> None:
    log: list[str] = []
    history = EditHistory(max_depth=2)
    for label in ("first", "second", "third"):
        history.push(_entry(label, log))
    assert history.undo() == "third"
    assert history.undo() == "second"
    with pytest.raises(IndexError, match="nothing to undo"):
        history.undo()


def test_clear_forgets_every_entry() -> None:
    log: list[str] = []
    history = EditHistory()
    history.push(_entry("first", log))
    history.push(_entry("second", log))
    history.undo()
    history.push(_entry("third", log))
    history.clear()
    assert not history.can_undo
    assert not history.can_redo
    with pytest.raises(IndexError, match="nothing to undo"):
        history.undo()


def test_empty_history_fails_closed() -> None:
    history = EditHistory()
    assert not history.can_undo
    assert not history.can_redo
    with pytest.raises(IndexError, match="nothing to undo"):
        history.undo()
    with pytest.raises(IndexError, match="nothing to redo"):
        history.redo()


def test_get_property_mirrors_the_routed_surface() -> None:
    application = Application(CFDXSession())
    application.set_property("numerics.cfl", 2.5)
    case = application.session.case
    assert get_property(case, "numerics.cfl") == 2.5
    assert get_property(case, "execution.solver") == (case.execution.solver or "")
    assert get_property(case, "physics.energy.enabled") is True
    with pytest.raises(KeyError, match="property is not editable"):
        get_property(case, "not.a.key")


def test_undo_restores_the_previous_cfl_value() -> None:
    application = Application(CFDXSession())
    initial = application.session.case.numerics.get("cfl", 1.0)
    application.set_property("numerics.cfl", 2.5)
    assert application.session.case.numerics["cfl"] == 2.5
    assert application.can_undo()
    assert not application.can_redo()
    description = application.undo()
    assert description == "set numerics.cfl to 2.5"
    assert application.session.case.numerics.get("cfl", 1.0) == initial
    assert not application.can_undo()
    assert application.can_redo()
    assert application.redo() == "set numerics.cfl to 2.5"
    assert application.session.case.numerics["cfl"] == 2.5


def test_successive_edits_undo_in_reverse_order() -> None:
    application = Application(CFDXSession())
    initial = application.session.case.numerics.get("cfl", 1.0)
    application.set_property("numerics.cfl", 2.5)
    application.set_property("numerics.cfl", 3.5)
    assert application.undo() == "set numerics.cfl to 3.5"
    assert application.session.case.numerics.get("cfl", 1.0) == 2.5
    assert application.undo() == "set numerics.cfl to 2.5"
    assert application.session.case.numerics.get("cfl", 1.0) == initial
    assert application.redo() == "set numerics.cfl to 2.5"
    assert application.session.case.numerics["cfl"] == 2.5


def test_application_edit_clears_the_redo_branch() -> None:
    application = Application(CFDXSession())
    application.set_property("numerics.cfl", 2.5)
    application.undo()
    assert application.can_redo()
    application.set_property("numerics.cfl", 3.5)
    assert not application.can_redo()
    assert application.can_undo()


def test_unchanged_value_is_not_recorded() -> None:
    application = Application(CFDXSession())
    current = application.session.case.numerics.get("cfl", 1.0)
    application.set_property("numerics.cfl", current)
    assert not application.can_undo()
    assert not application.can_redo()


def test_undo_restores_a_physics_toggle() -> None:
    application = Application(CFDXSession())
    application.session.case.physics["energy"] = {"enabled": True}
    application.set_property("physics.energy.enabled", False)
    assert application.session.case.physics["energy"]["enabled"] is False
    assert application.undo() == "set physics.energy.enabled to False"
    assert application.session.case.physics["energy"]["enabled"] is True


def test_new_project_clears_the_history() -> None:
    application = Application(CFDXSession())
    application.set_property("numerics.cfl", 2.5)
    assert application.can_undo()
    application.new_project()
    assert not application.can_undo()
    assert not application.can_redo()


def test_undo_publishes_the_application_state() -> None:
    application = Application(CFDXSession())
    published: list[ApplicationStateChanged] = []
    application.events.subscribe(ApplicationStateChanged, published.append)
    application.set_property("numerics.cfl", 2.5)
    application.undo()
    assert len(published) == 4
    assert published[-1].state.dirty
