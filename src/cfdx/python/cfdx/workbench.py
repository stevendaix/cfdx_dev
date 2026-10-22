"""Dockable CFDX Workbench composition root.

Domain actions remain owned by the application/session layer; Qt wires those
contracts to the optional mesh and result renderer adapters.
"""

from __future__ import annotations

from pathlib import Path

from .application import (
    Application,
    ApplicationStateChanged,
    ResultsChanged,
    WorkflowStatus,
    workflow_children,
)
from .mesh_model import read_mesh_catalog
from .session import CFDXSession

try:
    from PySide6.QtCore import QSettings, Qt, Signal
    from PySide6.QtGui import QAction, QBrush, QColor
    from PySide6.QtWidgets import (
        QDockWidget,
        QFileDialog,
        QFormLayout,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QMainWindow,
        QPushButton,
        QStackedLayout,
        QStatusBar,
        QToolBar,
        QTreeWidget,
        QTreeWidgetItem,
        QVBoxLayout,
        QWidget,
    )
except ImportError:  # pragma: no cover
    QMainWindow = object


if QMainWindow is not object:
    from .gui_3d import PyVistaQtView
    from .mesh_browser_panel import MeshBrowserPanel
    from .probes_panel import ProbesPanel
    from .results_panel import ResultsPanel
    from .run_center_panel import RunCenterPanel
    from .setup_panel import CaseSetupPanel

    class CFDXWorkbenchWindow(QMainWindow):
        """Initial Workbench composition root with stable dock object names."""

        # Execution callbacks originate from SolverRunner worker threads. Qt
        # widgets must only be touched on the GUI thread, so application events
        # cross this boundary through queued Qt signal delivery.
        _state_event = Signal(object)
        _results_event = Signal(object)

        SETTINGS_ORGANIZATION = "CFDX"
        SETTINGS_APPLICATION = "Workbench"

        def __init__(
            self,
            session: CFDXSession | None = None,
            *,
            application: Application | None = None,
        ) -> None:
            super().__init__()
            self.application = application or Application(session or CFDXSession())
            self.session = self.application.session
            self._application_state = self.application.state
            self._state_event.connect(self._state_changed)
            self._results_event.connect(self._results_changed)
            self.application.events.subscribe(
                ApplicationStateChanged, self._state_event.emit
            )
            self.application.events.subscribe(ResultsChanged, self._results_event.emit)
            self.setWindowTitle(f"CFDX Workbench — {self.session.case.name}")
            self.resize(1440, 900)
            self.setDockNestingEnabled(True)
            self._settings = QSettings(
                self.SETTINGS_ORGANIZATION, self.SETTINGS_APPLICATION
            )
            self._build_toolbar()
            self._build_docks()
            self._build_status_bar()
            self._restore_layout()

        def _build_toolbar(self) -> None:
            toolbar = QToolBar("Workbench", self)
            toolbar.setObjectName("workbench.toolbar")
            toolbar.setMovable(False)
            self.addToolBar(Qt.ToolBarArea.TopToolBarArea, toolbar)

            for label in (
                "New Project",
                "Open",
                "Save",
                "Open Mesh",
                "Open DAT Result",
                "Open DAT Checkpoint",
            ):
                action = QAction(label, self)
                action.setObjectName(
                    f"workbench.action.{label.lower().replace(' ', '_')}"
                )
                toolbar.addAction(action)
                if label == "New Project":
                    action.triggered.connect(self._new_project)
                    self.new_project_action = action
                elif label == "Open":
                    action.triggered.connect(self._open_case)
                    self.open_case_action = action
                elif label == "Open Mesh":
                    action.triggered.connect(self._open_mesh)
                    self.open_mesh_action = action
                elif label == "Open DAT Result":
                    action.triggered.connect(
                        lambda _checked=False: self._open_dat(use_for_restart=False)
                    )
                    self.open_dat_result_action = action
                elif label == "Open DAT Checkpoint":
                    action.triggered.connect(
                        lambda _checked=False: self._open_dat(use_for_restart=True)
                    )
                    self.open_dat_checkpoint_action = action
                else:
                    action.triggered.connect(self._save_case)
                    self.save_case_action = action
            toolbar.addSeparator()
            for label in ("Check", "Run", "Stop"):
                action = QAction(label, self)
                action.setObjectName(f"workbench.action.{label.lower()}")
                toolbar.addAction(action)
                if label == "Check":
                    action.triggered.connect(self._check_case)
                    self.check_action = action
                elif label == "Run":
                    action.triggered.connect(self._run)
                    self.run_action = action
                else:
                    action.triggered.connect(self._stop)
                    self.stop_action = action
            toolbar.addSeparator()
            for label in ("Contour", "Slice", "Glyph"):
                action = QAction(label, self)
                action.setObjectName(f"workbench.action.{label.lower()}")
                action.setEnabled(False)
                mode = label.lower()
                action.triggered.connect(
                    lambda _checked=False, mode=mode: self._postprocess(mode)
                )
                toolbar.addAction(action)
                setattr(self, f"{label.lower()}_action", action)
            toolbar.addSeparator()
            search = QAction("Search", self)
            search.setObjectName("workbench.action.search")
            search.setEnabled(False)
            toolbar.addAction(search)
            self.toolbar = toolbar

        def _new_project(self) -> None:
            self.application.new_project()
            self.session = self.application.session
            self._refresh_setup_case()
            self.setWindowTitle(f"CFDX Workbench — {self.session.case.name}")

        def _open_case(self) -> None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Open CFDX Case", "", "CFDX Case (*.cfdx.h5)"
            )
            if not path:
                return
            try:
                self.application.open_project(path)
                self.session = self.application.session
                self._refresh_setup_case()
                self.setWindowTitle(f"CFDX Workbench — {self.session.case.name}")
            except (OSError, ValueError) as exc:
                self.statusBar().showMessage(f"Open failed: {exc}")

        def _save_case(self) -> None:
            path = (
                Path(self.application.project_path)
                if self.application.project_path
                else None
            )
            if path is None:
                selected, _ = QFileDialog.getSaveFileName(
                    self,
                    "Save CFDX Case",
                    f"{self.session.case.name}.cfdx.h5",
                    "CFDX Case (*.cfdx.h5)",
                )
                if not selected:
                    return
                path = Path(selected)
            try:
                self.application.save_project(path)
                if self.application.restart_dat is not None:
                    self.statusBar().showMessage(
                        f"Case and restart checkpoint saved: {path.name}"
                    )
                else:
                    self.statusBar().showMessage(f"Case saved: {path.name}")
            except (OSError, ValueError) as exc:
                self.statusBar().showMessage(f"Save failed: {exc}")

        def _check_case(self) -> None:
            try:
                state = self.application.validate_case()
                if state.diagnostics:
                    self.statusBar().showMessage(
                        f"Validation: {len(state.diagnostics)} diagnostic(s)"
                    )
                else:
                    self.statusBar().showMessage("Validation: case is valid")
            except (OSError, TypeError, ValueError) as exc:
                self.statusBar().showMessage(f"Validation failed: {exc}")

        def _run(self) -> None:
            try:
                self.application.run()
            except (OSError, RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"Run failed: {exc}")

        def _stop(self) -> None:
            try:
                self.application.stop()
            except (RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"Stop failed: {exc}")

        def _build_docks(self) -> None:
            self.workflow_tree = QTreeWidget()
            self.workflow_tree.setObjectName("workbench.workflow_tree")
            self.workflow_tree.setHeaderLabel("Workflow")
            self.workflow_tree.itemSelectionChanged.connect(self._selection_changed)
            self._add_dock(
                "Workflow",
                "workbench.dock.workflow",
                self.workflow_tree,
                Qt.DockWidgetArea.LeftDockWidgetArea,
            )
            self._refresh_workflow(self._application_state)
            self.mesh_browser = MeshBrowserPanel()
            self.mesh_browser.selection_changed.connect(self._mesh_selection_changed)
            self._add_dock(
                "Mesh",
                "workbench.dock.mesh",
                self.mesh_browser,
                Qt.DockWidgetArea.LeftDockWidgetArea,
            )
            self.results_panel = ResultsPanel(self.application)
            self.results_panel.on_open = self._open_results_directory
            self._add_dock(
                "Results",
                "workbench.dock.results",
                self.results_panel,
                Qt.DockWidgetArea.LeftDockWidgetArea,
            )

            self.probes_panel = ProbesPanel(self.application)
            self.probes_panel.setObjectName("workbench.probes_panel")
            self._add_dock(
                "Probes",
                "workbench.dock.probes",
                self.probes_panel,
                Qt.DockWidgetArea.LeftDockWidgetArea,
            )

            self.setup_panel = CaseSetupPanel(self.session.case)
            self.setup_panel.setObjectName("workbench.setup_panel")
            self.setup_panel.changed.connect(self._setup_changed)
            self._add_dock(
                "Case Setup",
                "workbench.dock.setup",
                self.setup_panel,
                Qt.DockWidgetArea.RightDockWidgetArea,
            )

            try:
                renderer = PyVistaQtView()
                self.view3d = renderer
            except RuntimeError as exc:
                renderer = QLabel(f"3D renderer unavailable\n\n{exc}")
                self.view3d = None
            renderer.setObjectName("workbench.renderer")
            renderer.setMinimumSize(480, 320)

            self.viewport_placeholder = QWidget()
            self.viewport_placeholder.setObjectName("workbench.viewport.empty_state")
            self.viewport_placeholder.setStyleSheet("background: #f4f7fb;")
            empty_layout = QVBoxLayout(self.viewport_placeholder)
            empty_layout.setContentsMargins(32, 32, 32, 32)
            empty_layout.addStretch(1)
            empty_title_text = (
                "3D renderer unavailable"
                if self.view3d is None
                else "Your simulation workspace"
            )
            empty_title = QLabel(empty_title_text)
            empty_title.setObjectName("workbench.viewport.empty_title")
            empty_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
            empty_title.setStyleSheet(
                "font-size: 22px; font-weight: 600; color: #1f2937; "
                "background: transparent;"
            )
            empty_layout.addWidget(empty_title)
            empty_message_text = (
                str(renderer.text())
                if self.view3d is None
                else "Open a mesh or load a result to explore your CFD case in 3D."
            )
            empty_message = QLabel(empty_message_text)
            empty_message.setObjectName("workbench.viewport.empty_message")
            empty_message.setWordWrap(True)
            empty_message.setAlignment(Qt.AlignmentFlag.AlignCenter)
            empty_message.setStyleSheet(
                "font-size: 13px; color: #64748b; background: transparent;"
            )
            empty_layout.addWidget(empty_message)
            empty_actions = QWidget()
            empty_actions.setObjectName("workbench.viewport.empty_actions")
            actions_layout = QHBoxLayout(empty_actions)
            actions_layout.setContentsMargins(0, 12, 0, 0)
            actions_layout.setSpacing(12)
            for label, object_name, callback, primary in (
                ("Open Mesh", "open_mesh", self._open_mesh, True),
                (
                    "Load DAT Result",
                    "open_dat_result",
                    lambda: self._open_dat(use_for_restart=False),
                    False,
                ),
                (
                    "Open Results Folder",
                    "open_results",
                    self._open_results_directory,
                    False,
                ),
            ):
                button = QPushButton(label)
                button.setObjectName(f"workbench.empty_state.{object_name}")
                button.setMinimumHeight(38)
                button.setCursor(Qt.CursorShape.PointingHandCursor)
                button.setStyleSheet(
                    "QPushButton { padding: 8px 16px; border: 1px solid #cbd5e1; "
                    "border-radius: 6px; background: white; color: #334155; } "
                    "QPushButton:hover { background: #e2e8f0; } "
                    + (
                        "QPushButton { background: #2563eb; color: white; "
                        "border-color: #2563eb; font-weight: 600; } "
                        "QPushButton:hover { background: #1d4ed8; }"
                        if primary
                        else ""
                    )
                )
                button.clicked.connect(callback)
                actions_layout.addWidget(button)
            empty_layout.addWidget(empty_actions)
            empty_layout.addStretch(1)

            self.viewport_stack = QStackedLayout()
            self.viewport_stack.addWidget(self.viewport_placeholder)
            self.viewport_stack.addWidget(renderer)
            viewport = QWidget()
            viewport.setObjectName("workbench.viewport")
            viewport.setLayout(self.viewport_stack)
            self.setCentralWidget(viewport)

            properties = self._build_properties_panel()
            self._add_dock(
                "Properties",
                "workbench.dock.properties",
                properties,
                Qt.DockWidgetArea.RightDockWidgetArea,
            )

            monitor = RunCenterPanel(self.application)
            monitor.setObjectName("workbench.monitor")
            self._add_dock(
                "Monitors / Console",
                "workbench.dock.monitor",
                monitor,
                Qt.DockWidgetArea.BottomDockWidgetArea,
            )

        def _build_properties_panel(self) -> QWidget:
            panel = QWidget()
            panel.setObjectName("workbench.properties")
            panel._form = QFormLayout(panel)
            panel._form.addRow(QLabel("Select a workflow item"))
            self.properties_panel = panel
            self._refresh_properties()
            return panel

        def _open_results_directory(self) -> None:
            directory = QFileDialog.getExistingDirectory(
                self, "Open CFDX Results Directory"
            )
            if not directory:
                return
            try:
                self.application.open_results(directory)
            except (OSError, ValueError) as exc:
                self.statusBar().showMessage(f"Results error: {exc}")

        def _open_mesh(self) -> None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Open CFDX Mesh", "", "CFDX mesh (*.h5);;All files (*)"
            )
            if not path:
                return
            try:
                self.mesh_browser.set_catalog(read_mesh_catalog(Path(path)))
                if self.view3d is not None:
                    self.view3d.load_cfdx_mesh(path)
                    self._show_viewport()
                self.statusBar().showMessage(f"Mesh loaded: {path}")
            except (OSError, RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"Mesh error: {exc}")

        def _open_dat(self, *, use_for_restart: bool) -> None:
            if self.application.project_path is None:
                self.statusBar().showMessage(
                    "Save or open a CFDX case before loading a DAT file"
                )
                return
            path, _ = QFileDialog.getOpenFileName(
                self,
                "Open DAT Checkpoint" if use_for_restart else "Open DAT Result",
                str(Path(self.application.project_path).parent),
                "CFDX DAT (*.dat *.dat.h5 *.h5);;All files (*)",
            )
            if not path:
                return
            try:
                loaded = self.application.load_dat(
                    path, use_for_restart=use_for_restart
                )
                fields = []
                if self.view3d is not None:
                    fields = self.view3d.load_cfdx_dat(
                        str(self.application.project_path), path
                    )
                    self._show_viewport()
                self.statusBar().showMessage(
                    f"DAT loaded: {Path(path).name} | iteration={loaded.iteration} | "
                    f"time={loaded.time:g} | fields={len(fields or loaded.fields)}"
                    + (" | restart enabled" if use_for_restart else "")
                )
            except (OSError, RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"DAT error: {exc}")

        def _results_changed(self, event) -> None:
            result = event.state.results.selected_frame
            if self.view3d is None or result is None or not result.complete:
                return
            try:
                self.view3d.load(result.path)
                if event.state.results.selected_field:
                    self.view3d.set_field(event.state.results.selected_field)
                self._show_viewport()
            except (OSError, RuntimeError, ValueError, KeyError) as exc:
                self.statusBar().showMessage(f"Renderer error: {exc}")
            self._refresh_postprocess_actions(event.state)

        def _show_viewport(self) -> None:
            """Reveal the renderer after a dataset has loaded successfully."""
            if self.view3d is not None:
                self.viewport_stack.setCurrentWidget(self.view3d)

        def _refresh_postprocess_actions(self, state) -> None:
            enabled = (
                self.view3d is not None and state.results.selected_frame is not None
            )
            for name in ("contour_action", "slice_action", "glyph_action"):
                if hasattr(self, name):
                    getattr(self, name).setEnabled(enabled)

        def _postprocess(self, mode: str) -> None:
            if self.view3d is None:
                return
            field = self._application_state.results.selected_field
            try:
                if mode == "contour":
                    if not field:
                        raise ValueError("select a scalar result field before Contour")
                    self.view3d.contour(field)
                elif mode == "slice":
                    self.view3d.slice()
                elif mode == "glyph":
                    if not field:
                        raise ValueError("select a vector result field before Glyph")
                    self.view3d.glyph(field)
                self.statusBar().showMessage(f"Post-processing: {mode}")
            except (KeyError, RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"Post-processing error: {exc}")

        def _mesh_selection_changed(self, selection) -> None:
            if self.view3d is not None and selection.kind == "patch":
                try:
                    self.view3d.select(
                        selection.stable_id or f"patch:{selection.index}"
                    )
                except KeyError:
                    pass

        def _refresh_properties(self) -> None:
            if not hasattr(self, "properties_panel"):
                return
            form = self.properties_panel._form
            while form.count():
                item = form.takeAt(0)
                if item.widget() is not None:
                    item.widget().deleteLater()
            properties = self.application.properties()
            if not properties:
                form.addRow(QLabel("No editable properties for this selection"))
                return
            for prop in properties:
                editor = QLineEdit(str(prop.value))
                editor.setObjectName(f"workbench.property.{prop.key}")
                editor.setReadOnly(not prop.editable)
                button = QPushButton("Apply")
                button.setEnabled(prop.editable)
                button.clicked.connect(
                    lambda _checked=False, p=prop, e=editor: self._apply_property(
                        p.key, e.text()
                    )
                )
                form.addRow(prop.label, editor)
                form.addRow("", button)

        def _setup_changed(self) -> None:
            """Propagate specialized setup edits through the application state."""
            self.application.set_project_path(self.application.project_path, dirty=True)

        def _refresh_setup_case(self) -> None:
            if hasattr(self, "setup_panel"):
                self.setup_panel.set_case(self.session.case)

        def _apply_property(self, key: str, raw_value: str) -> None:
            try:
                value = float(raw_value) if key == "numerics.cfl" else raw_value
                self.application.set_property(key, value)
            except (KeyError, TypeError, ValueError) as exc:
                self.statusBar().showMessage(f"Property error: {exc}")
            else:
                self._refresh_properties()

        def _add_dock(
            self, title: str, object_name: str, widget: QWidget, area: Qt.DockWidgetArea
        ) -> QDockWidget:
            dock = QDockWidget(title, self)
            dock.setObjectName(object_name)
            dock.setWidget(widget)
            dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
            self.addDockWidget(area, dock)
            return dock

        def _build_status_bar(self) -> None:
            status = QStatusBar(self)
            status.setObjectName("workbench.status_bar")
            status.showMessage(
                f"State: {self.session.state.value} | Iteration: {self.session.iteration} | Time: {self.session.time:g}"
            )
            self.setStatusBar(status)

        @staticmethod
        def _status_color(status: WorkflowStatus) -> str:
            return {
                WorkflowStatus.COMPLETE: "#15803d",
                WorkflowStatus.WARNING: "#b45309",
                WorkflowStatus.ERROR: "#b91c1c",
                WorkflowStatus.NOT_CONFIGURED: "#64748b",
                WorkflowStatus.CONDITIONAL: "#7c3aed",
            }[status]

        @staticmethod
        def _status_marker(status: WorkflowStatus) -> str:
            return {
                WorkflowStatus.COMPLETE: "[OK]",
                WorkflowStatus.WARNING: "[!]",
                WorkflowStatus.ERROR: "[X]",
                WorkflowStatus.NOT_CONFIGURED: "[ ]",
                # Conditional means declared but not currently producible; never
                # render it as complete or silently omit the workflow capability.
                WorkflowStatus.CONDITIONAL: "[~]",
            }[status]

        def _refresh_workflow(self, state) -> None:
            self.workflow_tree.clear()
            roots = workflow_children(state.workflow, None)
            for step in roots:
                root = QTreeWidgetItem(
                    [f"{step.label}  {self._status_marker(step.status)}"]
                )
                root.setData(0, Qt.ItemDataRole.UserRole, step.id)
                root.setData(0, Qt.ItemDataRole.UserRole + 1, step.kind)
                root.setToolTip(0, step.message or step.status.value)
                root.setForeground(0, QBrush(QColor(self._status_color(step.status))))
                root_font = root.font(0)
                root_font.setBold(True)
                root.setFont(0, root_font)
                self.workflow_tree.addTopLevelItem(root)
                for child in workflow_children(state.workflow, step.id):
                    item = QTreeWidgetItem(
                        [f"{child.label}  {self._status_marker(child.status)}"]
                    )
                    item.setData(0, Qt.ItemDataRole.UserRole, child.id)
                    item.setData(0, Qt.ItemDataRole.UserRole + 1, child.kind)
                    item.setToolTip(0, child.message or child.status.value)
                    item.setForeground(
                        0, QBrush(QColor(self._status_color(child.status)))
                    )
                    root.addChild(item)
                root.setExpanded(True)

        def _state_changed(self, event) -> None:
            state = event.state
            self._application_state = state
            self._refresh_workflow(state)
            self.statusBar().showMessage(
                f"State: {state.simulation_state.value} | Iteration: {state.execution.iteration} | "
                f"Time: {state.execution.time:g}"
            )
            self._refresh_properties()
            if hasattr(self, "setup_panel"):
                self.setup_panel.set_diagnostics(state.diagnostics)
            if hasattr(self, "run_action"):
                running = state.simulation_state.value == "RUNNING"
                paused = state.simulation_state.value == "PAUSED"
                self.run_action.setEnabled(not running and not paused)
                self.stop_action.setEnabled(running or paused)
                self.check_action.setEnabled(not running and not paused)
                self.save_case_action.setEnabled(not running and not paused)

        def _selection_changed(self) -> None:
            items = self.workflow_tree.selectedItems()
            if not items:
                return
            item = items[0]
            stable_id = item.data(0, Qt.ItemDataRole.UserRole)
            kind = item.data(0, Qt.ItemDataRole.UserRole + 1)
            label = item.text(0).rsplit("  [", 1)[0]
            self.application.select(stable_id, kind, label)

        def _restore_layout(self) -> None:
            geometry = self._settings.value("geometry")
            state = self._settings.value("dock_state")
            if geometry is not None:
                self.restoreGeometry(geometry)
            if state is not None:
                self.restoreState(state)

        def closeEvent(self, event) -> None:
            self._settings.setValue("geometry", self.saveGeometry())
            self._settings.setValue("dock_state", self.saveState())
            super().closeEvent(event)

else:

    class CFDXWorkbenchWindow:
        """Clear failure when optional GUI dependencies are unavailable."""

        def __init__(self, *args, **kwargs) -> None:
            raise RuntimeError("PySide6 is required for the CFDX GUI")
