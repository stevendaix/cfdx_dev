"""Stable capability registry for the CFDX Scientific Workbench.

Gap G3 of the Phase 0 audit (docs/gui/PHASE0_AUDIT.md): actions must never be
silently disabled. Every capability exposes a stable identifier matching its
Qt action objectName, a human label, its roadmap lot, and prerequisites that
are evaluated against a plain context object. This module has no Qt import so
availability rules stay testable headlessly and reusable by CLI/TUI clients.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Prerequisite:
    key: str
    description: str


PREREQUISITES = {
    "renderer": Prerequisite("renderer", "a 3D renderer"),
    "dataset": Prerequisite("dataset", "a mesh or result loaded in the 3D viewport"),
    "idle": Prerequisite("idle", "no simulation currently running or paused"),
    "active": Prerequisite("active", "a running or paused simulation"),
}


@dataclass(frozen=True)
class CapabilityContext:
    """Facts availability is computed from; never holds domain state."""

    renderer_ready: bool = False
    dataset_loaded: bool = False
    simulation_active: bool = False


@dataclass(frozen=True)
class Capability:
    id: str
    label: str
    lot: int
    prerequisites: tuple[Prerequisite, ...]

    @property
    def prerequisites_tip(self) -> str:
        if not self.prerequisites:
            return ""
        names = " and ".join(p.description for p in self.prerequisites)
        return f"Requires: {names}"


@dataclass(frozen=True)
class CapabilityStatus:
    capability: Capability
    available: bool
    missing: tuple[Prerequisite, ...]

    @property
    def reason(self) -> str:
        if not self.missing:
            return ""
        names = " and ".join(p.description for p in self.missing)
        return f"requires {names}"

    @property
    def tip(self) -> str:
        """Clear-text tip naming only the currently missing prerequisites."""
        if self.available:
            return ""
        names = " and ".join(p.description for p in self.missing)
        return f"Requires: {names}"


_EVALUATORS = {
    "renderer": lambda context: context.renderer_ready,
    "dataset": lambda context: context.dataset_loaded,
    "idle": lambda context: not context.simulation_active,
    "active": lambda context: context.simulation_active,
}


CAPABILITIES: dict[str, Capability] = {
    capability.id: capability
    for capability in (
        Capability(
            "workbench.action.contour",
            "Contour",
            3,
            (PREREQUISITES["renderer"], PREREQUISITES["dataset"]),
        ),
        Capability(
            "workbench.action.slice",
            "Slice",
            3,
            (PREREQUISITES["renderer"], PREREQUISITES["dataset"]),
        ),
        Capability(
            "workbench.action.glyph",
            "Glyph",
            3,
            (PREREQUISITES["renderer"], PREREQUISITES["dataset"]),
        ),
        Capability(
            "workbench.action.run",
            "Run",
            1,
            (PREREQUISITES["idle"],),
        ),
        Capability(
            "workbench.action.check",
            "Check",
            1,
            (PREREQUISITES["idle"],),
        ),
        Capability(
            "workbench.action.save",
            "Save",
            1,
            (PREREQUISITES["idle"],),
        ),
        Capability(
            "workbench.action.stop",
            "Stop",
            1,
            (PREREQUISITES["active"],),
        ),
    )
}


def get_capability(capability_id: str) -> Capability | None:
    return CAPABILITIES.get(capability_id)


def evaluate_capability(
    capability_id: str, context: CapabilityContext
) -> CapabilityStatus:
    capability = CAPABILITIES.get(capability_id)
    if capability is None:
        raise KeyError(f"unknown capability: {capability_id}")
    missing = tuple(
        prerequisite
        for prerequisite in capability.prerequisites
        if not _EVALUATORS[prerequisite.key](context)
    )
    return CapabilityStatus(capability, not missing, missing)
