"""Explicit application commands used by GUI, TUI, and CLI frontends."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol, TYPE_CHECKING

from ..session import CFDXSession, ChangeImpact

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
