"""Application facade shared by CFDX frontends."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from ..dat_io import DatRestart, read_dat_restart
from ..execution import ExecutionController
from ..mesh_model import MeshCatalog, read_mesh_catalog
from ..probe import ProbeCatalog
from ..probe_validation import ProbeSeries, read_probe_csv, write_probe_csv
from ..project import Project
from ..results_series import ResultSeries, discover_result_series
from ..runner import SolverRunner
from ..session import CFDXSession, ChangeImpact
from ..setup_schema import SetupSchema, build_setup_schema
from ..validation import validate_case as _validate_case
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
from .properties import PropertyState, properties_for_selection
from .results import ResultsState, build_results_state
from .run_center import RunCenterModel
from .state import ApplicationState, SelectionState, build_application_state


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
        self._restart_dat: Path | None = None
        self.mesh_catalog: MeshCatalog | None = None
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

    def attach_controller(self, controller: ExecutionController) -> ApplicationState:
        if self.run_center is not None:
            self.run_center.close()
        self.controller = controller
        self.run_center = RunCenterModel(controller)
        self.run_center.on_change = lambda _state: self._publish_state()
        return self._publish_state()

    def replace_session(
        self, session: CFDXSession, *, project_path: str | Path | None = None
    ) -> ApplicationState:
        runner = (
            getattr(self.controller, "runner", None)
            if self.controller is not None
            else None
        )
        if getattr(runner, "running", False):
            raise RuntimeError("cannot replace an active session")
        if self.run_center is not None:
            self.run_center.close()
        self.session = session
        self.controller = None
        self.run_center = None
        self.project_path = project_path
        self.selection = SelectionState()
        self.diagnostics = ()
        self.results = ResultsState()
        self._restart_dat = None
        self.mesh_catalog = None
        self.dirty = False
        return self._publish_state()

    def open_project(
        self,
        path: str | Path,
        *,
        with_dat: bool = False,
        dat_path: str | Path | None = None,
    ) -> ApplicationState:
        """Open a canonical CFDX project through the shared case lifecycle."""
        project = Project(Path(path))
        restart_path: Path | None = None
        if with_dat:
            session, restart_path = project.load_with_dat(
                Path(dat_path) if dat_path is not None else None
            )
        else:
            session = project.load()
        state = self.replace_session(session, project_path=project.path)
        if restart_path is not None:
            self._restart_dat = restart_path
            state = self._publish_state()
        return state

    @property
    def restart_dat(self) -> Path | None:
        """Return the explicitly selected numerical checkpoint, if any."""
        return self._restart_dat

    def new_project(self) -> ApplicationState:
        """Create a new in-memory project through the shared lifecycle."""
        return self.replace_session(CFDXSession())

    def save_project(self, path: str | Path | None = None) -> ApplicationState:
        """Save the current setup/mesh state without consuming DAT state."""
        target = Path(path) if path is not None else self.project_path
        if target is None:
            raise ValueError("a .cfdx.h5 project path is required")
        project = Project(Path(target))
        if self._restart_dat is not None:
            _, canonical_dat = project.save_with_dat(self.session, self._restart_dat)
            self._restart_dat = canonical_dat
        else:
            project.save(self.session)
        return self.set_project_path(project.path, dirty=False)

    def save_project_with_dat(
        self, source_dat: str | Path, path: str | Path | None = None
    ) -> ApplicationState:
        """Save setup plus an explicit numerical checkpoint as a paired DAT."""
        target = Path(path) if path is not None else self.project_path
        if target is None:
            raise ValueError("a .cfdx.h5 project path is required")
        project = Project(Path(target))
        _, canonical_dat = project.save_with_dat(self.session, Path(source_dat))
        self._restart_dat = canonical_dat
        return self.set_project_path(project.path, dirty=False)

    def load_dat(
        self, path: str | Path, *, use_for_restart: bool = False
    ) -> DatRestart:
        """Read a numerical DAT artifact without involving a frontend renderer."""
        if self.project_path is None:
            raise ValueError("load a CFDX case before loading a DAT file")
        restart = read_dat_restart(Path(path))
        if use_for_restart:
            self._restart_dat = Path(path)
            self._publish_state()
        return restart

    def clear_restart(self) -> ApplicationState:
        """Clear the selected numerical restart without changing case setup."""
        self._restart_dat = None
        return self._publish_state()

    def load_restart(
        self,
        dat_path: str | Path,
        *,
        solver: str | None = None,
        mesh: str | Path | None = None,
        output_dir: str | Path | None = None,
    ) -> DatRestart:
        """Load a bare DAT checkpoint as the solver restart state.

        The production solver writes a standalone ``restart.dat`` next to its
        output directory; the GUI must be able to open that file directly and
        have it become the solver's restart state, without requiring the
        caller to first save a paired ``.cfdx.h5`` project. The DAT is read
        through the existing ``dat_io`` contract, so the numerical values are
        validated exactly as they would be for a paired checkpoint.

        The production solver requires ``--mesh`` and ``--output-dir`` even
        when restarting from a checkpoint, so a bare DAT restart must carry
        those options alongside the checkpoint path. They are stored on the
        session's execution config rather than on the DAT itself, because the
        DAT carries numerical state only.
        """
        path = Path(dat_path)
        if not path.is_file():
            raise FileNotFoundError(path)
        restart = read_dat_restart(path)
        if solver is not None:
            self.session.case.execution.solver = solver
        elif not self.session.case.execution.solver:
            raise ValueError(
                "a bare DAT restart requires a configured execution.solver"
            )
        if mesh is not None:
            self.session.case.execution.mesh = str(Path(mesh))
        if output_dir is not None:
            self.session.case.execution.output_dir = str(Path(output_dir))
        self._restart_dat = path
        self._publish_state()
        return restart

    def ensure_controller(self) -> ExecutionController:
        """Construct the canonical solver execution controller on demand.

        Two paths are supported:

        * A saved ``.cfdx.h5`` project: the solver command carries the
          project path and, if a paired DAT was loaded, the restart option and
          checkpoint path.
        * A bare DAT restart loaded through :meth:`load_restart`: the solver
          command carries the DAT path directly through the case-configured
          restart option. The DAT has no case setup of its own, so the command
          is built from the DAT path and the configured solver rather than
          from a project path.
        """
        if self.controller is not None:
            return self.controller
        solver = self.session.case.execution.solver
        if not solver:
            raise ValueError("execution.solver must be configured before Run")
        restart_option = self.session.case.execution.restart_option
        if self._restart_dat is not None and not restart_option:
            raise ValueError(
                "a DAT checkpoint is loaded but no restart option is configured"
            )
        if self.project_path is not None:
            command: list[str] = [solver, str(self.project_path)]
            cwd = Path(self.project_path).parent
            if self._restart_dat is not None:
                command.extend([restart_option, str(self._restart_dat)])
        elif self._restart_dat is not None:
            # A bare DAT restart: the solver consumes the checkpoint
            # directly, so the command carries the DAT path through the
            # restart option alongside the mesh and output-dir the solver
            # requires even when restarting from a checkpoint.
            command = [solver, "--mesh", str(self.session.case.execution.mesh)]
            command += ["--output-dir", str(self.session.case.execution.output_dir)]
            command += [restart_option, str(self._restart_dat)]
            cwd = Path(self._restart_dat).parent
        else:
            raise ValueError("save the case before starting the solver")
        probe_catalog = ProbeCatalog(self.session.case.probes)
        if probe_catalog.probes:
            for spec in probe_catalog.specs():
                command += ["--probe", spec]
            command += [
                "--probe-csv",
                str(probe_catalog.csv_path(self.project_path or self._restart_dat)),
            ]
        if self.session.case.execution.mpi_ranks > 1:
            command = [
                "mpiexec",
                "-n",
                str(self.session.case.execution.mpi_ranks),
                *command,
            ]
        controller = ExecutionController(
            self.session,
            SolverRunner(command, cwd=cwd),
        )
        self.attach_controller(controller)
        return controller

    def validate_case(self, mesh=None) -> ApplicationState:
        """Validate the current case through the shared application contract."""
        return self.validate(_validate_case, mesh)

    def set_project_path(
        self, path: str | Path | None, *, dirty: bool | None = None
    ) -> ApplicationState:
        self.project_path = path
        if dirty is not None:
            self.dirty = dirty
        return self._publish_state()

    def set_numerical_option(
        self, key: str, value: Any, impact: ChangeImpact = ChangeImpact.HOT
    ) -> ApplicationState:
        return self.execute(SetNumericalOption(key, value, impact))

    def setup_schema(self) -> SetupSchema:
        """Return the headless schema for the current case."""
        return build_setup_schema(self.session.case)

    def properties(self) -> tuple[PropertyState, ...]:
        return properties_for_selection(
            self.session.case,
            self.selection.stable_id,
            catalog=self.mesh_catalog,
        )

    def set_property(self, key: str, value: Any) -> ApplicationState:
        return self.execute(SetProperty(key, value))

    def set_diagnostics(self, diagnostics: Iterable[Any]) -> ApplicationState:
        self.diagnostics = tuple(diagnostics)
        return self._publish_state()

    def validate(self, validator, mesh=None) -> ApplicationState:
        report = validator(self.session.case, mesh)
        return self.set_diagnostics(report.diagnostics)

    def read_probe_results(self, path: str | Path) -> tuple[ProbeSeries, ...]:
        """Read native solver probe history through the shared application API."""
        return read_probe_csv(Path(path))

    def export_probe_results(
        self, path: str | Path, series: Iterable[ProbeSeries]
    ) -> Path:
        """Export probe history using the canonical native CSV schema."""
        return write_probe_csv(Path(path), tuple(series))

    def set_results(
        self, series: ResultSeries | None, *, directory: str | Path | None = None
    ) -> ApplicationState:
        self.results = build_results_state(series, directory=directory)
        self.events.publish(ResultsChanged(self.state))
        return self.state

    def open_results(self, directory: str | Path) -> ApplicationState:
        path = Path(directory)
        return self.set_results(
            discover_result_series(path, inspect_fields=True), directory=path
        )

    def load_mesh(self, path: str | Path) -> MeshCatalog:
        """Load a CFDX mesh and keep its catalog for contextual properties."""
        self.mesh_catalog = read_mesh_catalog(Path(path))
        return self.mesh_catalog

    def read_dat_restart(self, path: str | Path) -> DatRestart:
        """Read a numerical DAT artifact via the application facade."""
        return read_dat_restart(Path(path))

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
        """Start the configured solver through the shared application facade."""
        self.ensure_controller()
        return self.execute(RunSolver())

    def pause(self) -> ApplicationState:
        return self.execute(PauseSolver())

    def resume(self) -> ApplicationState:
        return self.execute(ResumeSolver())

    def stop(self) -> ApplicationState:
        return self.execute(StopSolver())

    def select(
        self, stable_id: str | None, kind: str | None, label: str | None
    ) -> ApplicationState:
        self.selection = SelectionState(stable_id, kind, label)
        snapshot = self.state
        self.events.publish(SelectionChanged(snapshot))
        return snapshot
