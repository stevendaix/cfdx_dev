"""Explicit application commands used by GUI, TUI, and CLI frontends."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol, TYPE_CHECKING

from ..probe import ProbeCatalog
from ..session import CFDXSession, ChangeImpact
from .properties import set_property

if TYPE_CHECKING:
    from ..execution import ExecutionController


class CommandContext(Protocol):
    session: CFDXSession
    controller: ExecutionController | None


class Command:
    """A user intent that is independent from Qt widgets."""

    def execute(self, context: CommandContext) -> None:
        raise NotImplementedError


@dataclass(frozen=True)
class SetNumericalOption(Command):
    key: str
    value: Any
    impact: ChangeImpact = ChangeImpact.HOT

    def execute(self, context: CommandContext) -> None:
        context.session.edit(self.key, self.value, self.impact)


@dataclass(frozen=True)
class SetProperty(Command):
    key: str
    value: Any

    def execute(self, context: CommandContext) -> None:
        impact = set_property(context.session.case, self.key, self.value)
        context.session.case_revision += 1
        if impact is not ChangeImpact.HOT:
            context.session.mark_restart_required()


@dataclass(frozen=True)
class ConfigureProbes(Command):
    """Replace the case's probe catalogue as a reversible setup action.

    Adding or removing probes only reconfigures solver instrumentation (it does
    not alter the discretisation), so the impact is ``HOT`` and no restart is
    required. The catalogue is forwarded to the solver as ``--probe``/``--probe``
    arguments by the execution backend.
    """

    catalog: ProbeCatalog

    def execute(self, context: CommandContext) -> None:
        self.catalog.validate()
        context.session.case.probes = list(self.catalog.probes)
        context.session.case_revision += 1


class _ControllerCommand(Command):
    def _controller(self, context: CommandContext) -> ExecutionController:
        if context.controller is None:
            raise RuntimeError("an execution controller is required")
        return context.controller


class RunSolver(_ControllerCommand):
    def execute(self, context: CommandContext) -> None:
        self._controller(context).start()


class PauseSolver(_ControllerCommand):
    def execute(self, context: CommandContext) -> None:
        self._controller(context).pause()


class ResumeSolver(_ControllerCommand):
    def execute(self, context: CommandContext) -> None:
        self._controller(context).resume()


class StopSolver(_ControllerCommand):
    def execute(self, context: CommandContext) -> None:
        self._controller(context).stop()


class ValidateCase(Command):
    def __init__(self, validator, mesh=None) -> None:
        self.validator = validator
        self.mesh = mesh
        self.report = None

    def execute(self, context: CommandContext) -> None:
        self.report = self.validator(context.session.case, self.mesh)
