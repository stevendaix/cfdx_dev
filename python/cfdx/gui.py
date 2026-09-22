"""Optional PySide6 GUI for CFDX case setup and solver execution."""
from __future__ import annotations

from pathlib import Path

from .case_io import read_case, read_case_with_dat, save_case, save_case_with_dat
from .execution import ExecutionController
from .output import OutputThrottle
from .runner import SolverRunner
from .session import CFDXSession, SimulationState, ChangeImpact
from .tui import TuiRenderer
from .watcher import ResultWatcher
from .mesh_model import read_mesh_catalog
from .mesh_browser_panel import MeshBrowserPanel
from .results_series import discover_result_series
from .results_series_panel import ResultsSeriesPanel

try:
    from PySide6.QtCore import QObject, QTimer, Qt, Signal
    from PySide6.QtWidgets import (
        QApplication, QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout,
        QLabel, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton,
        QTabWidget, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
    )
except ImportError:  # pragma: no cover
    QApplication = None
    QObject = object
    Signal = None


if QApplication is not None:
    from .setup_panel import CaseSetupPanel
    from .gui_3d import PyVistaQtView

    class SessionSignals(QObject):
        state_changed = Signal(object)
        output = Signal(str, bool)
        metrics_changed = Signal(object)
        results_changed = Signal()

    class CFDXMainWindow(QMainWindow):
        """Fluent-like shell with one session, controller and explicit case files."""

        def __init__(
            self,
            session: CFDXSession | None = None,
            results_dir: Path | None = None,
            execution_controller: ExecutionController | None = None,
        ) -> None:
            super().__init__()
            self.session = session or CFDXSession()
            self.controller = execution_controller
            self._case_path: Path | None = None
            self._restart_dat: Path | None = None
            self._result_source: Path | None = None
            self._dirty = False
            self.signals = SessionSignals()

            self.setWindowTitle(f"CFDX — {self.session.case.name}")
            self.resize(1200, 760)

            self._output_throttle = OutputThrottle(self._append_output_now)
            self._flush_timer = QTimer(self)
            self._flush_timer.setInterval(50)
            self._flush_timer.timeout.connect(self._output_throttle.flush)
            self._flush_timer.start()
            self._result_watcher = (
                ResultWatcher(results_dir, self.signals.results_changed.emit)
                if results_dir is not None else None
            )

            self._build_actions()
            self._build_ui()
            self._connect_controller()
            self.refresh()
            if self._result_watcher is not None:
                self._result_watcher.start()

        def _build_actions(self) -> None:
            file_menu = self.menuBar().addMenu("File")
            self.open_case_action = file_menu.addAction("Read Case…")
            self.open_case_dat_action = file_menu.addAction("Read Case + DAT…")
            file_menu.addSeparator()
            self.save_case_action = file_menu.addAction("Save Case")
            self.save_case_as_action = file_menu.addAction("Save Case As…")
            self.save_case_dat_action = file_menu.addAction("Save Case + DAT…")
            file_menu.addSeparator()
            self.open_mesh_action = file_menu.addAction("Open Mesh…")
            self.open_results_action = file_menu.addAction("Open Results Directory…")
            file_menu.addSeparator()
            self.exit_action = file_menu.addAction("Exit")

            run_menu = self.menuBar().addMenu("Run")
            self.run_action = run_menu.addAction("Run")
            self.pause_action = run_menu.addAction("Pause")
            self.resume_action = run_menu.addAction("Resume")
            self.stop_action = run_menu.addAction("Stop")

            self.open_case_action.triggered.connect(self._open_case)
            self.open_case_dat_action.triggered.connect(self._open_case_with_dat)
            self.save_case_action.triggered.connect(self._save_case)
            self.save_case_as_action.triggered.connect(self._save_case_as)
            self.save_case_dat_action.triggered.connect(self._save_case_with_dat)
            self.open_mesh_action.triggered.connect(self._open_mesh)
            self.open_results_action.triggered.connect(self._open_results_directory)
            self.exit_action.triggered.connect(self.close)
            self.run_action.triggered.connect(self._run)
            self.pause_action.triggered.connect(self._pause)
            self.resume_action.triggered.connect(self._resume)
            self.stop_action.triggered.connect(self._stop)

        def _build_ui(self) -> None:
            central = QWidget()
            layout = QHBoxLayout(central)
            self.tree = QTreeWidget()
            self.tree.setHeaderLabels(["Case"])
            self.log = QPlainTextEdit()
            self.log.setReadOnly(True)
            self.status = QLabel()
            self.file_status = QLabel("Unsaved case")
            self.file_status.setToolTip("Persistent case and restart artifacts")

            self.run_button = QPushButton("Run")
            self.pause_button = QPushButton("Pause")
            self.resume_button = QPushButton("Resume")
            self.stop_button = QPushButton("Stop")
            controls = QHBoxLayout()
            controls.addWidget(self.run_button)
            controls.addWidget(self.pause_button)
            controls.addWidget(self.resume_button)
            controls.addWidget(self.stop_button)

            self.parameters = QFormLayout()
            self.cfl = QDoubleSpinBox()
            self.cfl.setRange(0.0, 1.0e6)
            self.cfl.setDecimals(6)
            self.cfl.setValue(float(self.session.case.numerics.get("cfl", 1.0)))
            self.cfl.valueChanged.connect(self._set_cfl)
            self.parameters.addRow("CFL", self.cfl)

            self.mesh_browser = MeshBrowserPanel()
            self.mesh_browser.selection_changed.connect(self._mesh_selection_changed)
            self.results_series = ResultsSeriesPanel()
            self.results_series.frame_changed.connect(self._result_frame_changed)

            self.setup_panel = CaseSetupPanel(self.session.case)
            self.setup_panel.changed.connect(self._mark_dirty)
            right = QVBoxLayout()
            right.addWidget(self.status)
            right.addWidget(self.file_status)
            right.addLayout(controls)
            right.addLayout(self.parameters)
            right.addWidget(self.mesh_browser)
            right.addWidget(self.setup_panel)
            right.addWidget(self.log)

            setup_container = QWidget()
            setup_container.setLayout(right)

            visualization = QWidget()
            visualization_layout = QVBoxLayout(visualization)
            result_controls = QHBoxLayout()
            self.open_result_button = QPushButton("Open VTU / VTK…")
            self.result_status = QLabel("No result loaded")
            result_controls.addWidget(self.open_result_button)
            result_controls.addWidget(self.result_status, 1)
            visualization_layout.addLayout(result_controls)
            visualization_layout.addWidget(self.results_series)
            try:
                self.view3d = PyVistaQtView()
                visualization_layout.addWidget(self.view3d, 1)
            except RuntimeError as exc:
                self.view3d = None
                self.result_status.setText(str(exc))

            tabs = QTabWidget()
            tabs.addTab(setup_container, "Case / Run")
            tabs.addTab(visualization, "3D Results")

            layout.addWidget(self.tree, 1)
            layout.addWidget(tabs, 3)
            self.setCentralWidget(central)
            self.open_result_button.clicked.connect(self._open_result)

            self.run_button.clicked.connect(self._run)
            self.pause_button.clicked.connect(self._pause)
            self.resume_button.clicked.connect(self._resume)
            self.stop_button.clicked.connect(self._stop)
            self.signals.state_changed.connect(self._refresh_status)
            self.signals.output.connect(self._queue_output)
            self.signals.metrics_changed.connect(self._refresh_metrics)
            self.signals.results_changed.connect(self.refresh)

        def _open_mesh(self) -> bool:
            path, _ = QFileDialog.getOpenFileName(
                self, "Open CFDX Mesh", str(self._case_path.parent if self._case_path else ""),
                "CFDX mesh (*.h5);;All files (*)"
            )
            if not path:
                return False
            try:
                catalog = read_mesh_catalog(Path(path))
                self.mesh_browser.set_catalog(catalog)
                if self.view3d is not None:
                    self.view3d.load_cfdx_mesh(path)
                return True
            except (OSError, ValueError) as exc:
                self._show_error("Open Mesh failed", str(exc))
                return False

        def _open_results_directory(self) -> bool:
            directory = QFileDialog.getExistingDirectory(
                self, "Open CFDX Results Directory",
                str(self._case_path.parent if self._case_path else "")
            )
            if not directory:
                return False
            try:
                series = discover_result_series(Path(directory))
                self.results_series.set_series(series)
                if series.frames:
                    self._result_frame_changed(series.frames[0])
                else:
                    self.result_status.setText("No VTK result files found")
                return True
            except (OSError, ValueError) as exc:
                self._show_error("Open Results failed", str(exc))
                return False

        def _result_frame_changed(self, frame) -> None:
            if self.view3d is None:
                self.result_status.setText(str(frame.path))
                return
            try:
                self.view3d.load(str(frame.path))
                self._result_source = frame.path
                self.result_status.setText(str(frame.path))
            except (OSError, RuntimeError, ValueError) as exc:
                self._show_error("Open 3D result failed", str(exc))

        def _mesh_selection_changed(self, selection) -> None:
            if self.view3d is not None and selection.kind == "patch":
                try:
                    self.view3d.select(f"patch:{selection.index}")
                except KeyError:
                    pass
            self.result_status.setText(
                f"Selected {selection.kind} {selection.index}"
                + (f" ({selection.name})" if selection.name else "")
            )

        def _open_result(self) -> bool:
            if self.view3d is None:
                self._show_error("3D view unavailable", self.result_status.text())
                return False
            path, _ = QFileDialog.getOpenFileName(
                self,
                "Open CFDX Result",
                str(self._case_path.parent if self._case_path else ""),
                "VTK results (*.vtu *.vtk *.vtp *.pvtu);;All files (*)",
            )
            if not path:
                return False
            try:
                self.view3d.load(path)
                self._result_source = Path(path)
                self.result_status.setText(str(self._result_source))
                return True
            except (OSError, RuntimeError, ValueError) as exc:
                self._show_error("Open 3D result failed", str(exc))
                return False

        def _connect_controller(self) -> None:
            if self.controller is None:
                return
            self.controller.on_output = lambda line, stderr: self.signals.output.emit(line, stderr)
            self.controller.on_metrics = lambda metrics: self.signals.metrics_changed.emit(metrics)
            self.controller.on_complete = lambda _result: self.signals.state_changed.emit(self.session.state)

        def _solver_command(self) -> list[str]:
            solver = self.session.case.execution.solver
            if not solver:
                raise ValueError("execution.solver must be configured before Run")
            if self._case_path is None:
                raise ValueError("Save the case before Run")
            command = [solver, str(self._case_path)]
            if self._restart_dat is not None:
                restart_option = self.session.case.execution.restart_option
                if not restart_option:
                    raise ValueError("A DAT restart is loaded but no solver restart option is configured")
                command.extend([restart_option, str(self._restart_dat)])
            if self.session.case.execution.mpi_ranks > 1:
                command = [
                    "mpiexec", "-n", str(self.session.case.execution.mpi_ranks), *command                ]
            return command

        def _ensure_controller(self) -> ExecutionController:
            if self.controller is None:
                self.controller = ExecutionController(
                    self.session, SolverRunner(self._solver_command(), cwd=self._case_path.parent)
                )
                self._connect_controller()
            return self.controller

        def _ensure_case_path(self) -> bool:
            if self._case_path is not None:
                return True
            return self._save_case_as()

        def _mark_dirty(self) -> None:
            self._dirty = True
            self._refresh_file_status()

        def _save_case_as(self) -> bool:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save CFDX Case", self.session.case.name + ".cfdx.h5",
                "CFDX Case (*.cfdx.h5 *.h5)"
            )
            if not path:
                return False
            try:
                self._case_path = save_case(self.session, Path(path))
                self._restart_dat = None
                self._dirty = False
                self._refresh_file_status()
                return True
            except (OSError, ValueError) as exc:
                self._show_error("Save Case failed", str(exc))
                return False

        def _save_case(self) -> bool:
            if self._case_path is None:
                return self._save_case_as()
            try:
                save_case(self.session, self._case_path)
                self._dirty = False
                self._refresh_file_status()
                return True
            except (OSError, ValueError) as exc:
                self._show_error("Save Case failed", str(exc))
                return False

        def _save_case_with_dat(self) -> bool:
            if not self._ensure_case_path():
                return False
            dat_path, _ = QFileDialog.getOpenFileName(
                self, "Select solver DAT restart", str(self._case_path.parent),
                "Solver restart (*.dat);;All files (*)"
            )
            if not dat_path:
                return False
            try:
                self._case_path, self._restart_dat = save_case_with_dat(
                    self.session, self._case_path, Path(dat_path)
                )
                self._dirty = False
                self._refresh_file_status()
                return True
            except (OSError, ValueError) as exc:
                self._show_error("Save Case + DAT failed", str(exc))
                return False

        def _open_case(self) -> bool:
            path, _ = QFileDialog.getOpenFileName(
                self, "Read CFDX Case", "", "CFDX Case (*.cfdx.h5 *.h5)"
            )
            if not path:
                return False
            try:
                loaded = read_case(Path(path))
                self._replace_session(loaded)
                self._case_path = Path(path)
                self._restart_dat = None
                self._dirty = False
                self._refresh_file_status()
                return True
            except (OSError, ValueError) as exc:
                self._show_error("Read Case failed", str(exc))
                return False

        def _open_case_with_dat(self) -> bool:
            path, _ = QFileDialog.getOpenFileName(
                self, "Read CFDX Case + DAT", "", "CFDX Case (*.cfdx.h5 *.h5)"
            )
            if not path:
                return False
            dat_path, _ = QFileDialog.getOpenFileName(
                self, "Select paired solver DAT", str(Path(path).parent),
                "Solver restart (*.dat);;All files (*)"
            )
            if not dat_path:
                return False
            try:
                loaded, restart = read_case_with_dat(Path(path), Path(dat_path))
                self._replace_session(loaded)
                self._case_path = Path(path)
                self._restart_dat = restart
                self._dirty = False
                self._refresh_file_status()
                return True
            except (OSError, ValueError) as exc:
                self._show_error("Read Case + DAT failed", str(exc))
                return False

        def _replace_session(self, session: CFDXSession) -> None:
            if self.controller is not None and self.controller.runner.running:
                raise RuntimeError("Cannot replace an active session")
            self.session = session
            self.setWindowTitle(f"CFDX — {self.session.case.name}")
            self.setup_panel.set_case(self.session.case)
            try:
                self.setup_panel.changed.disconnect(self._mark_dirty)
            except (RuntimeError, TypeError):
                pass
            self.setup_panel.changed.connect(self._mark_dirty)
            self.cfl.blockSignals(True)
            self.cfl.setValue(float(self.session.case.numerics.get("cfl", 1.0)))
            self.cfl.blockSignals(False)
            self.controller = None
            self._connect_controller()
            self.refresh()

        def refresh(self) -> None:
            self.tree.clear()
            for node in self.session.case_tree():
                item = QTreeWidgetItem([node.label])
                item.setData(0, Qt.ItemDataRole.UserRole, node.id)
                self.tree.addTopLevelItem(item)
            self._refresh_file_status()
            self._refresh_status(self.session.state)

        def _refresh_metrics(self, _metrics) -> None:
            self._refresh_status(self.session.state)

        def _refresh_file_status(self) -> None:
            if self._case_path is None:
                self.file_status.setText("Unsaved case" + (" • Modified" if self._dirty else ""))
                return
            restart = f" | DAT: {self._restart_dat.name}" if self._restart_dat else ""
            self.file_status.setText(f"Case: {self._case_path}{restart}" + (" • Modified" if self._dirty else ""))

        def _refresh_status(self, state: SimulationState, *_args) -> None:
            self.status.setText(
                f"State: {state.value} | Iteration: {self.session.iteration} | "
                f"Time: {self.session.time:g}"
            )
            run_enabled = state in {
                SimulationState.CREATED, SimulationState.READY,
                SimulationState.PAUSED, SimulationState.STOPPED,
            }
            self.run_button.setEnabled(run_enabled)
            self.run_action.setEnabled(run_enabled)
            self.pause_button.setEnabled(state is SimulationState.RUNNING)
            self.pause_action.setEnabled(state is SimulationState.RUNNING)
            self.resume_button.setEnabled(state is SimulationState.PAUSED)
            self.resume_action.setEnabled(state is SimulationState.PAUSED)
            self.stop_button.setEnabled(state in {SimulationState.RUNNING, SimulationState.PAUSED})
            self.stop_action.setEnabled(state in {SimulationState.RUNNING, SimulationState.PAUSED})
            self.save_case_action.setEnabled(state not in {SimulationState.RUNNING, SimulationState.VALIDATING})
            self.save_case_as_action.setEnabled(state not in {SimulationState.RUNNING, SimulationState.VALIDATING})
            self.save_case_dat_action.setEnabled(state not in {SimulationState.RUNNING, SimulationState.VALIDATING})

        def _set_cfl(self, value: float) -> None:
            if self.session.state in {SimulationState.RUNNING, SimulationState.VALIDATING}:
                return
            self.session.edit("cfl", value, ChangeImpact.HOT)
            self._mark_dirty()

        def _run(self) -> None:
            try:
                if self._dirty:
                    choice = QMessageBox.question(
                        self,
                        "Unsaved changes",
                        "Save the case before starting the solver?",
                        QMessageBox.StandardButton.Save
                        | QMessageBox.StandardButton.Discard
                        | QMessageBox.StandardButton.Cancel,
                        QMessageBox.StandardButton.Save,
                    )
                    if choice is QMessageBox.StandardButton.Cancel:
                        return
                    if choice is QMessageBox.StandardButton.Save and not self._save_case():
                        return
                if not self._ensure_case_path():
                    return
                controller = self._ensure_controller()
                controller.start()
                self.signals.state_changed.emit(self.session.state)
            except (OSError, RuntimeError, ValueError, NotImplementedError) as exc:
                self._show_error("Run failed", str(exc))
                self.signals.state_changed.emit(self.session.state)

        def _pause(self) -> None:
            try:
                self._ensure_controller().pause()
                self.signals.state_changed.emit(self.session.state)
            except (RuntimeError, NotImplementedError) as exc:
                self._show_error("Pause failed", str(exc))

        def _resume(self) -> None:
            try:
                self._ensure_controller().resume()
                self.signals.state_changed.emit(self.session.state)
            except (RuntimeError, NotImplementedError) as exc:
                self._show_error("Resume failed", str(exc))

        def _stop(self) -> None:
            try:
                self._ensure_controller().stop()
                self.signals.state_changed.emit(self.session.state)
            except RuntimeError as exc:
                self._show_error("Stop failed", str(exc))

        def _show_error(self, title: str, message: str) -> None:
            QMessageBox.critical(self, title, message)

        def _queue_output(self, line: str, is_stderr: bool = False) -> None:
            self._output_throttle.push(line, is_stderr)

        def append_output(self, line: str, is_stderr: bool = False) -> None:
            self._queue_output(line, is_stderr)

        def _append_output_now(self, line: str, is_stderr: bool) -> None:
            self.log.appendPlainText(("[stderr] " if is_stderr else "") + line)

        def render_tui(self) -> str:
            return TuiRenderer.render(self.session)

        def closeEvent(self, event) -> None:
            if self._dirty and (self.controller is None or not self.controller.runner.running):
                choice = QMessageBox.question(
                    self,
                    "Unsaved changes",
                    "Save changes before closing?",
                    QMessageBox.StandardButton.Save
                    | QMessageBox.StandardButton.Discard
                    | QMessageBox.StandardButton.Cancel,
                    QMessageBox.StandardButton.Save,
                )
                if choice is QMessageBox.StandardButton.Cancel:
                    event.ignore()
                    return
                if choice is QMessageBox.StandardButton.Save and not self._save_case():
                    event.ignore()
                    return
            if self.controller is not None and self.controller.runner.running:
                choice = QMessageBox.question(
                    self, "Solver still running",
                    "The solver is still running. Stop it before closing?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.Yes,
                )
                if choice is QMessageBox.StandardButton.No:
                    event.ignore()
                    return
                self._stop()
            self._flush_timer.stop()
            if self._result_watcher is not None:
                self._result_watcher.stop()
            if self.view3d is not None:
                self.view3d.close()
            super().closeEvent(event)

    def create_application(argv: list[str] | None = None) -> QApplication:
        return QApplication.instance() or QApplication(argv or [])

    def launch(session: CFDXSession | None = None, argv: list[str] | None = None) -> int:
        app = create_application(argv)
        window = CFDXMainWindow(session)
        window.show()
        return app.exec()
else: