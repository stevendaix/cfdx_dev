"""Dependency-free TUI renderer for CFDX sessions."""
from __future__ import annotations

from .session import CFDXSession


class TuiRenderer:
    """Render application state without requiring a terminal UI framework."""

    @staticmethod
    def render(session: CFDXSession) -> str:
        lines = [
            f"CFDX | Case: {session.case.name}",
            f"State: {session.state.value}",
            f"Iteration: {session.iteration} | Time: {session.time:g}",
            f"Revision: case={session.case_revision} mesh={session.mesh_revision} "
            f"physics={session.physics_revision} numerics={session.numerics_revision}",
            "[RUN] [PAUSE] [STOP] [EDIT] [CHECKPOINT]",
        ]
        if session.requires_restart:
            lines.append("WARNING: restart/rebuild required before continuing.")
        return "\n".join(lines)
