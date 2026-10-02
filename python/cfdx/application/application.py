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
from .events import ApplicationStateChanged, EventBus, ResultsChanged, SelectionChanged
from .state import ApplicationState, SelectionState, build_application_state
from .properties import PropertyState, properties_for_selection
from .run_center import RunCenterModel
from .results import ResultsState, build_results_state
from ..results_series import ResultSeries, discover_result_series


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
        self.results = ResultsState()
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
            results=self.results,
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

    def set_results(self, series: ResultSeries | None, *, directory: str | Path | None = None) -> ApplicationState:
        self.results = build_results_state(series, directory=directory)
        self.events.publish(ResultsChanged(self.state))
        return self.state

    def open_results(self, directory: str | Path) -> ApplicationState:
        path = Path(directory)
        return self.set_results(discover_result_series(path, inspect_fields=True), directory=path)

    def select_result_frame(self, stable_id: str) -> ApplicationState:
        if stable_id not in {frame.stable_id for frame in self.results.frames}:
            raise KeyError(stable_id)
        self.results = ResultsState(
            self.results.directory,
            self.results.frames,
            stable_id,
            self.results.selected_field,
        )
        self.events.publish(ResultsChanged(self.state))
        return self.state

    def select_result_field(self, field: str) -> ApplicationState:
        if field not in self.results.field_names:
            raise KeyError(field)
        self.results = ResultsState(
            self.results.directory,
            self.results.frames,
            self.results.selected_frame_id,
            field,
        )
        self.events.publish(ResultsChanged(self.state))
        return self.state

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
