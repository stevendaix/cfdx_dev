"""Bounded undo/redo history of reversible application edits.

Gap G4 of the Phase 0 audit (docs/gui/PHASE0_AUDIT.md): the Workbench had no
undo. Reversible edits are recorded by the application layer itself — the
same facade serves GUI, CLI, and TUI — so every frontend inherits the same
history semantics. An entry carries the inverse and replay callables of one
edit; recording a new edit clears the redo branch, and the undo depth stays
bounded so a long session cannot grow without limit.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

MAX_UNDO_DEPTH = 100


@dataclass(frozen=True)
class HistoryEntry:
    """One reversible edit: a clear-text label plus inverse and replay."""

    description: str
    undo: Callable[[], None]
    redo: Callable[[], None]


class EditHistory:
    """Two stacks over reversible edits; a new edit clears the redo branch."""

    def __init__(self, max_depth: int = MAX_UNDO_DEPTH) -> None:
        self._max_depth = max_depth
        self._undo_entries: list[HistoryEntry] = []
        self._redo_entries: list[HistoryEntry] = []

    @property
    def can_undo(self) -> bool:
        return bool(self._undo_entries)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo_entries)

    def push(self, entry: HistoryEntry) -> None:
        """Record one reversible edit and drop any obsolete redo branch."""
        self._redo_entries.clear()
        self._undo_entries.append(entry)
        if len(self._undo_entries) > self._max_depth:
            self._undo_entries.pop(0)

    def undo(self) -> str:
        """Revert the most recent edit and return its clear-text label."""
        if not self._undo_entries:
            raise IndexError("nothing to undo")
        entry = self._undo_entries.pop()
        entry.undo()
        self._redo_entries.append(entry)
        return entry.description

    def redo(self) -> str:
        """Replay the most recent undone edit and return its label."""
        if not self._redo_entries:
            raise IndexError("nothing to redo")
        entry = self._redo_entries.pop()
        entry.redo()
        self._undo_entries.append(entry)
        return entry.description

    def clear(self) -> None:
        """Forget every reversible edit, e.g. when the session is replaced."""
        self._undo_entries.clear()
        self._redo_entries.clear()
