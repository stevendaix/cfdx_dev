"""Presentation mapping for solver and simulation states (P0 review).

The GUI must distinguish RUNNING, CONVERGED, STAGNATED, DIVERGED and STOPPED
at a glance, and RUNNING must never be presented as established convergence.
Source strings come from the session ``SimulationState`` or from the solver's
``to_string(convergence_status)``; the GUI never guesses. Unrecognized values
are surfaced verbatim through ``unknown_state_view`` instead of being coerced
to a nearby state. Token names reference theme TOKENS keys, not raw colors,
so this module stays Qt-free and testable headlessly.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SolverStateView:
    key: str
    label: str
    marker: str
    token: str
    converged: bool = False
    terminal: bool = False


SESSION_STATE_VIEWS: dict[str, SolverStateView] = {
    "CREATED": SolverStateView("CREATED", "Created", "○", "text.muted"),
    "VALIDATING": SolverStateView("VALIDATING", "Validating", "◌", "status.info"),
    "READY": SolverStateView("READY", "Ready", "✓", "status.success"),
    "RUNNING": SolverStateView(
        "RUNNING",
        "Running — convergence not established",
        "●",
        "accent.primary",
    ),
    "PAUSED": SolverStateView("PAUSED", "Paused", "⏸", "status.warning"),
    "STOPPED": SolverStateView(
        "STOPPED", "Stopped by user", "■", "status.warning", terminal=True
    ),
    "CONVERGED": SolverStateView(
        "CONVERGED",
        "Converged",
        "✓",
        "status.success",
        converged=True,
        terminal=True,
    ),
    "FAILED": SolverStateView("FAILED", "Failed", "✕", "status.error", terminal=True),
}

CONVERGENCE_STATUS_VIEWS: dict[str, SolverStateView] = {
    "RUNNING": SESSION_STATE_VIEWS["RUNNING"],
    "CONVERGED": SESSION_STATE_VIEWS["CONVERGED"],
    "STAGNATED": SolverStateView(
        "STAGNATED", "Stagnated", "≈", "status.warning", terminal=True
    ),
    "DIVERGED": SolverStateView(
        "DIVERGED", "Diverged", "↗", "status.error", terminal=True
    ),
    "MAX_ITERATIONS": SolverStateView(
        "MAX_ITERATIONS",
        "Maximum iterations reached",
        "#",
        "status.warning",
        terminal=True,
    ),
    "STOPPED": SESSION_STATE_VIEWS["STOPPED"],
}


def session_state_view(value: str) -> SolverStateView:
    view = SESSION_STATE_VIEWS.get(value)
    if view is None:
        return unknown_state_view(value)
    return view


def convergence_status_view(value: str) -> SolverStateView:
    view = CONVERGENCE_STATUS_VIEWS.get(value)
    if view is None:
        return unknown_state_view(value)
    return view


def unknown_state_view(value: str) -> SolverStateView:
    """Surface an unrecognized solver value verbatim instead of guessing."""
    text = str(value)
    return SolverStateView(text, text, "?", "status.warning")
