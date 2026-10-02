"""Application facade shared by CFDX frontends."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from ..execution import ExecutionController
from ..session import CFDXSession, ChangeImpact
from .commands import (
    Command,
    PauseSolver,
    ResumeSolver,
    RunSolver,
    SetNumericalOption,
    SetProperty,
    StopSolver,
)
from .events import ApplicationStateChanged, EventBus, SelectionChanged
from .state import ApplicationState, SelectionState, build_application_state
from .properties import PropertyState, properties_for_selection
from .run_center import RunCenterModel


class Application:
    """Application API boundary between frontends and CFDX orchestration."""

    def __init__(
        self,
        session: CFDXSession | None = None,
        *,
        controller: ExecutionController | None = None,
        project_path: str | Path | None = None,
        capabilities: Iterable[str] = (),
    ) -> None:
        self.session = session or CFDXSession()
        self.controller = controller
        self.project_path = project_path
        self.capabilities = tuple(capabilities)
        self.selection = SelectionState()
        self.diagnostics: tuple[Any, ...] = ()
        self.dirty = False
        self.events = EventBus()
        self.run_center = RunCenterModel(controller) if controller is not None else None
        if self.run_center is not None:
            self.run_center.on_change = lambda _state: self._publish_state()

    @property
    def state(self) -> ApplicationState:
        return build_application_state(
            self.session,
            project_path=self.project_path,
            dirty=self.dirty,
            selection=self.selection,
            controller=self.controller,
            capabilities=self.capabilities,
            diagnostics=self.diagnostics,
        )

    def execute(self, command: Command) -> ApplicationState:
        command.execute(self)
        if isinstance(command, (SetNumericalOption, SetProperty)):
            self.dirty = True
        return self._publish_state()

    def _publish_state(self) -> ApplicationState:
        snapshot = self.state
        self.events.publish(ApplicationStateChanged(snapshot))
        return snapshot

    def set_numerical_option(
        self, key: str, value: Any, impact: ChangeImpact = ChangeImpact.HOT
    ) -> ApplicationState:
        return self.execute(SetNumericalOption(key, value, impact))

    def properties(self) -> tuple[PropertyState, ...]:
        return properties_for_selection(self.session.case, self.selection.stable_id)

    def set_property(self, key: str, value: Any) -> ApplicationState:
        return self.execute(SetProperty(key, value))

    def set_diagnostics(self, diagnostics: Iterable[Any]) -> ApplicationState:
        self.diagnostics = tuple(diagnostics)
        return self._publish_state()

    def validate(self, validator, mesh=None) -> ApplicationState:
        report = validator(self.session.case, mesh)
        return self.set_diagnostics(report.diagnostics)

    def run(self) -> ApplicationState:
        return self.execute(RunSolver())

    def pause(self) -> ApplicationState:
        return self.execute(PauseSolver())

    def resume(self) -> ApplicationState:
        return self.execute(ResumeSolver())

    def stop(self) -> ApplicationState:
        return self.execute(StopSolver())

    def select(self, stable_id: str | None, kind: str | None, label: str | None) -> ApplicationState:
        self.selection = SelectionState(stable_id, kind, label)
        snapshot = self.state
        self.events.publish(SelectionChanged(snapshot))
        return snapshot
