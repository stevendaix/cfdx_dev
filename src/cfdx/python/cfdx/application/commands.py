"""Explicit application commands used by GUI, TUI, and CLI frontends."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Protocol

from ..probe import ProbeCatalog
from ..session import CFDXSession, ChangeImpact
from .properties import set_property

if TYPE_CHECKING:
    from ..execution import ExecutionController


class CommandContext(Protocol):
    session: CFDXSession
    controller: ExecutionController | None

    def ensure_controller(self) -> ExecutionController: ...


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
        if context.controller is not None:
            return context.controller
        # The Application owns controller construction through
        # ``ensure_controller``, which builds the canonical solver command
        # from the case config and, when a DAT restart is loaded, from the
        # checkpoint path. A controller command must therefore be able to
        # construct the controller on demand rather than requiring the
        # caller to have pre-built one.
        ensure = getattr(context, "ensure_controller", None)
        if callable(ensure):
            return ensure()
        raise RuntimeError("an execution controller is required")


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
