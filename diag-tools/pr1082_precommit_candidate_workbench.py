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
    SelectionChanged,
    WorkflowStatus,
    workflow_children,
)
from .capabilities import (
    Capability,
    CapabilityContext,
    evaluate_capability,
    get_capability,
)
from .selection import parse_selection_id, selection_label
from .session import CFDXSession
from .solver_states import session_state_view
from .theme import TOKENS, apply_workbench_theme

try:
    from PySide6.QtCore import QSettings, Qt, Signal
    from PySide6.QtGui import QAction, QActionGroup, QBrush, QColor
    from PySide6.QtWidgets import (
        QDialog,
        QDockWidget,
        QFileDialog,
        QHBoxLayout,
        QInputDialog,
        QLabel,
        QLineEdit,
        QListWidget,
        QListWidgetItem,
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
    from .properties_panel import PropertiesPanel
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
        _selection_event = Signal(object)

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
            self._selection_event.connect(self._selection_state_changed)
            self.application.events.subscribe(
                ApplicationStateChanged, self._state_event.emit
            )
            self.application.events.subscribe(ResultsChanged, self._results_event.emit)
            self.application.events.subscribe(
                SelectionChanged, self._selection_event.emit
            )
            self.setWindowTitle(f"CFDX Workbench — {self.session.case.name}")
            self.resize(1440, 900)
            self.setDockNestingEnabled(True)
            self._settings = QSettings(
                self.SETTINGS_ORGANIZATION, self.SETTINGS_APPLICATION
            )
            self._docks: dict[str, QDockWidget] = {}
            self._workspace = "setup"
            self._build_toolbar()
            self._build_docks()
            self._build_status_bar()
            apply_workbench_theme(self)
            self._restore_layout()
            self._set_workspace(str(self._settings.value("workspace", "setup")))
            self._refresh_action_availability()

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
                capability = get_capability(action.objectName())
                if capability is not None:
                    action.setStatusTip(capability.prerequisites_tip)
                mode = label.lower()
                action.triggered.connect(
                    lambda _checked=False, mode=mode: self._postprocess(mode)
                )
                toolbar.addAction(action)
                setattr(self, f"{label.lower()}_action", action)
            toolbar.addSeparator()
            search = QAction("Search", self)
            search.setObjectName("workbench.action.search")
            search.setShortcut("Ctrl+K")
            search.setToolTip("Search commands and workspaces (Ctrl+K)")
            search.triggered.connect(self._show_command_palette)
            toolbar.addAction(search)
            self.search_action = search
            self.toolbar = toolbar
            self._build_workspace_toolbar()
            self._register_capability_actions()

        def _build_workspace_toolbar(self) -> None:
            """Expose task-oriented layouts without duplicating domain actions."""
            self.addToolBarBreak(Qt.ToolBarArea.TopToolBarArea)
            toolbar = QToolBar("Workspace", self)
            toolbar.setObjectName("workbench.workspace_toolbar")
            toolbar.setMovable(False)
            group = QActionGroup(self)
            group.setExclusive(True)
            self.workspace_actions: dict[str, QAction] = {}
            for workspace, label in (
                ("setup", "Setup"),
                ("mesh", "Mesh"),
                ("run", "Run"),
                ("analyze", "Analyze"),
            ):
                action = QAction(label, self)
                action.setObjectName(f"workbench.workspace.{workspace}")
                action.setCheckable(True)
                group.addAction(action)
                toolbar.addAction(action)
                action.triggered.connect(
                    lambda _checked=False, key=workspace: self._set_workspace(key)
                )
                self.workspace_actions[workspace] = action
            self.addToolBar(Qt.ToolBarArea.TopToolBarArea, toolbar)
            self.workspace_toolbar = toolbar

        def _register_capability_actions(self) -> None:
            """Bind every registered capability id to its Qt action.

            The stable capability id equals the Qt action objectName, so the
            binding stays declarative: an action outside the registry never
            enters the palette, and tests guard the bijection.
            """
            self._capability_actions = {
                action.objectName(): action
                for action in self.findChildren(QAction)
                if get_capability(action.objectName()) is not None
            }

        def _show_command_palette(self) -> None:
            """Open a keyboard-first palette over the capability registry."""
            dialog = QDialog(self)
            dialog.setObjectName("workbench.command_palette")
            dialog.setWindowTitle("Command Palette")
            dialog.setMinimumWidth(520)
            layout = QVBoxLayout(dialog)
            query = QLineEdit(dialog)
            query.setObjectName("workbench.command_palette.query")
            query.setPlaceholderText("Search commands and workspaces...")
            results = QListWidget(dialog)
            results.setObjectName("workbench.command_palette.results")
            layout.addWidget(query)
            layout.addWidget(results)

            commands: list[tuple[Capability, QAction]] = []
            for capability_id, action in self._capability_actions.items():
                capability = get_capability(capability_id)
                if capability is not None:
                    commands.append((capability, action))
            commands.sort(key=lambda item: item[0].label.casefold())

            def refresh_results(text: str = "") -> None:
                results.clear()
                needle = text.casefold().strip()
                context = self._capability_context()
                for capability, action in commands:
                    if needle and needle not in capability.label.casefold():
                        continue
                    status = evaluate_capability(capability.id, context)
                    if status.available:
                        item = QListWidgetItem(capability.label)
                    else:
                        # P0 (Phase 0 audit, G3): never hide a disabled
                        # command silently; show why it is unavailable.
                        reason = status.tip or "not available in the current state"
                        item = QListWidgetItem(f"{capability.label} — {reason}")
                        item.setToolTip(reason)
                    item.setData(Qt.ItemDataRole.UserRole, action)
                    results.addItem(item)
                if results.count():
                    results.setCurrentRow(0)

            def activate_current() -> None:
                item = results.currentItem()
                if item is None:
                    return
                action = item.data(Qt.ItemDataRole.UserRole)
                dialog.accept()
                if isinstance(action, QAction) and action.isEnabled():
                    action.trigger()

            query.textChanged.connect(refresh_results)
            results.itemActivated.connect(lambda _item: activate_current())
            query.returnPressed.connect(activate_current)
            refresh_results()
            dialog.show()
            query.setFocus()
            self.command_palette_dialog = dialog
            self.command_palette_query = query
            self.command_palette_results = results

        def _set_workspace(self, workspace: str) -> None:
            """Apply a named dock layout while preserving all panel instances."""
            layouts = {
                "setup": {
                    "workbench.dock.workflow",
                    "workbench.dock.mesh",
                    "workbench.dock.setup",
                    "workbench.dock.properties",
                },
                "mesh": {
                    "workbench.dock.workflow",
                    "workbench.dock.mesh",
                    "workbench.dock.properties",
                },
                "run": {
                    "workbench.dock.workflow",
                    "workbench.dock.setup",
                    "workbench.dock.properties",
                    "workbench.dock.monitor",
                },
                "analyze": {
                    "workbench.dock.results",
                    "workbench.dock.probes",
                    "workbench.dock.visualization_pipeline",
                    "workbench.dock.properties",
                    "workbench.dock.monitor",
                },
            }
            if workspace not in layouts:
                workspace = "setup"
            visible_docks = layouts[workspace]
            for object_name, dock in self._docks.items():
                dock.setVisible(object_name in visible_docks)
            self._workspace = workspace
            for key, action in self.workspace_actions.items():
                action.setChecked(key == workspace)

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

            self.visualization_operations = QListWidget()
            self.visualization_operations.setObjectName(
                "workbench.visualization.operations"
            )
            self.remove_visualization_operation_button = QPushButton("Remove selected")
            self.remove_visualization_operation_button.setObjectName(
                "workbench.visualization.remove"
            )
            self.remove_visualization_operation_button.setEnabled(False)
            self.remove_visualization_operation_button.clicked.connect(
                self._remove_visualization_operation
            )
            self.reset_visualization_button = QPushButton("Reset pipeline")
            self.reset_visualization_button.setObjectName(
                "workbench.visualization.reset"
            )
            self.reset_visualization_button.setEnabled(False)
            self.reset_visualization_button.clicked.connect(
                self._reset_visualization_pipeline
            )
            visualization_panel = QWidget()
            visualization_layout = QVBoxLayout(visualization_panel)
            visualization_layout.addWidget(self.visualization_operations)
            visualization_layout.addWidget(self.remove_visualization_operation_button)
            visualization_layout.addWidget(self.reset_visualization_button)
            self.visualization_operations.currentRowChanged.connect(
                self._visualization_operation_selection_changed
            )
            self._add_dock(
                "Visualization Pipeline",
                "workbench.dock.visualization_pipeline",
                visualization_panel,
                Qt.DockWidgetArea.RightDockWidgetArea,
            )
            self._refresh_visualization_operations()

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
            title_font = empty_title.font()
            title_font.setPointSize(18)
            title_font.setBold(True)
            empty_title.setFont(title_font)
            empty_title.setAutoFillBackground(False)
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
            empty_message.setProperty("role", "secondary")
            empty_message.setAutoFillBackground(False)
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
                button.setProperty("emphasis", "primary" if primary else "secondary")
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
            """Compose the typed contextual inspector over the selection."""
            self.properties_panel = PropertiesPanel()
            self.properties_panel.on_apply = self._apply_property
            return self.properties_panel

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
                catalog = self.application.load_mesh(path)
                self.mesh_browser.set_catalog(catalog)
                if self.view3d is not None:
                    self.view3d.load_cfdx_mesh(path)
                    self._show_viewport()
                    self._refresh_action_availability()
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
                    self._refresh_action_availability()
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
            self._application_state = event.state
            self._refresh_action_availability()
            self._refresh_visualization_operations()

        def _show_viewport(self) -> None:
            """Reveal the renderer after a dataset has loaded successfully."""
            if self.view3d is not None:
                self.viewport_stack.setCurrentWidget(self.view3d)

        def _capability_context(self) -> CapabilityContext:
            """Snapshot the facts capability availability is computed from."""
            state = self._application_state
            view3d = getattr(self, "view3d", None)
            has_renderer = view3d is not None
            has_frame = state.results.selected_frame is not None
            if view3d is None:
                has_dataset = False
            else:
                has_dataset = (
                    view3d._dat_dataset is not None
                    or view3d._current_result is not None
                )
            simulation_active = state.simulation_state.value in ("RUNNING", "PAUSED")
            return CapabilityContext(
                renderer_ready=has_renderer,
                dataset_loaded=has_frame or has_dataset,
                simulation_active=simulation_active,
            )

        def _refresh_action_availability(self) -> None:
            """Compute every registered action's availability from the registry."""
            context = self._capability_context()
            for capability_id, action in self._capability_actions.items():
                status = evaluate_capability(capability_id, context)
                action.setEnabled(status.available)
                action.setStatusTip(status.tip)

        def _postprocess_fields(self, *, vector: bool) -> list[str]:
            """Return fields compatible with the requested visualization operation."""
            if self.view3d is None:
                return []
            dataset = self.view3d._dataset()
            names: set[str] = set()
            for attributes in (dataset.point_data, dataset.cell_data):
                for name in attributes:
                    values = attributes[name]
                    shape = getattr(values, "shape", ())
                    is_vector = len(shape) == 2 and shape[1] == 3
                    if is_vector == vector:
                        names.add(str(name))
            return sorted(names, key=str.casefold)

        def _postprocess(self, mode: str) -> None:
            """Collect operation-specific inputs before applying a visualization filter."""
            if getattr(self, "view3d", None) is None:
                return
            try:
                if mode in {"contour", "glyph"}:
                    fields = self._postprocess_fields(vector=mode == "glyph")
                    if not fields:
                        kind = "vector" if mode == "glyph" else "scalar"
                        self.statusBar().showMessage(
                            f"No compatible {kind} fields are available for {mode.title()}"
                        )
                        return
                    selected = self._application_state.results.selected_field
                    default_index = fields.index(selected) if selected in fields else 0
                    field, accepted = QInputDialog.getItem(
                        self,
                        f"{mode.title()} field",
                        "Vector field" if mode == "glyph" else "Scalar field",
                        fields,
                        default_index,
                        False,
                    )
                    if not accepted:
                        return
                    if mode == "contour":
                        isosurfaces, accepted = QInputDialog.getInt(
                            self,
                            "Contour settings",
                            "Number of isosurfaces:",
                            10,
                            2,
                            100,
                        )
                        if not accepted:
                            return
                        self.view3d.contour(field, isosurfaces=isosurfaces)
                    else:
                        factor, accepted = QInputDialog.getDouble(
                            self,
                            "Glyph settings",
                            "Scale factor:",
                            1.0,
                            0.001,
                            1000.0,
                            3,
                        )
                        if not accepted:
                            return
                        self.view3d.glyph(field, factor=factor)
                elif mode == "slice":
                    axes = ["X", "Y", "Z"]
                    axis, accepted = QInputDialog.getItem(
                        self, "Slice orientation", "Normal axis:", axes, 0, False
                    )
                    if not accepted:
                        return
                    normals = {
                        "X": (1.0, 0.0, 0.0),
                        "Y": (0.0, 1.0, 0.0),
                        "Z": (0.0, 0.0, 1.0),
                    }
                    self.view3d.slice(normal=normals[axis])
                else:
                    raise ValueError(f"unsupported visualization operation: {mode}")
                self.statusBar().showMessage(f"Post-processing: {mode}")
                self._refresh_visualization_operations()
            except (KeyError, RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"Post-processing error: {exc}")

        def _visualization_operation_selection_changed(self, row: int) -> None:
            self.remove_visualization_operation_button.setEnabled(
                row >= 0 and getattr(self, "view3d", None) is not None
            )

        def _refresh_visualization_operations(self) -> None:
            """Reflect the renderer pipeline in the Workbench operation list."""
            if not hasattr(self, "visualization_operations"):
                return
            self.visualization_operations.clear()
            view3d = getattr(self, "view3d", None)
            operations = (
                view3d._visualization_pipeline.operations if view3d is not None else ()
            )
            for index, operation in enumerate(operations, start=1):
                details = operation.kind
                if operation.kind == "contour":
                    details += (
                        f" — {operation.scalars}, {operation.isosurfaces} isosurfaces"
                    )
                elif operation.kind == "glyph":
                    details += f" — {operation.scalars}, scale {operation.factor:g}"
                elif operation.kind == "slice":
                    details += f" — normal {operation.normal}"
                item = QListWidgetItem(f"{index}. {details}")
                item.setData(Qt.ItemDataRole.UserRole, index - 1)
                self.visualization_operations.addItem(item)
            self.remove_visualization_operation_button.setEnabled(False)
            self.reset_visualization_button.setEnabled(bool(operations))

        def _remove_visualization_operation(self) -> None:
            """Remove the selected operation and replay the remaining pipeline."""
            if self.view3d is None:
                return
            item = self.visualization_operations.currentItem()
            if item is None:
                return
            index = item.data(Qt.ItemDataRole.UserRole)
            try:
                self.view3d.remove_visualization_operation(index)
            except (IndexError, RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"Visualization update failed: {exc}")
                return
            self._refresh_visualization_operations()
            self.statusBar().showMessage("Visualization operation removed")

        def _reset_visualization_pipeline(self) -> None:
            """Clear all derived operations and return to the original dataset."""
            if self.view3d is None:
                return
            try:
                self.view3d.reset_visualization()
            except (RuntimeError, ValueError) as exc:
                self.statusBar().showMessage(f"Visualization reset failed: {exc}")
                return
            self._refresh_visualization_operations()
            self.statusBar().showMessage("Visualization pipeline reset")

        def _mesh_selection_changed(self, selection) -> None:
            """Sync the 3D view, the inspector and the selection state."""
            if selection.kind not in {"patch", "face", "region"}:
                return
            stable_id = selection.stable_id
            if stable_id is None and selection.kind == "patch":
                stable_id = f"patch:{selection.index}"
            if stable_id is None:
                return
            try:
                summary = selection_label(parse_selection_id(stable_id))
            except ValueError:
                return
            self.application.select(stable_id, selection.kind, selection.name)
            if self.view3d is not None:
                try:
                    self.view3d.select(stable_id)
                    summary = self.view3d.selection_summary()
                except (KeyError, ValueError):
                    pass
            self.statusBar().showMessage(f"Selected: {summary}")

        def _refresh_properties(self) -> None:
            """Rebuild the typed inspector for the current selection context."""
            if not hasattr(self, "properties_panel"):
                return
            try:
                properties = self.application.properties()
            except (KeyError, ValueError):
                properties = ()
            self.properties_panel.set_properties(properties)

        def _setup_changed(self) -> None:
            """Propagate specialized setup edits through the application state."""
            self.application.set_project_path(self.application.project_path, dirty=True)

        def _refresh_setup_case(self) -> None:
            if hasattr(self, "setup_panel"):
                self.setup_panel.set_case(self.session.case)

        def _apply_property(self, key: str, value: object) -> None:
            """Apply one typed property value and report failures in clear text."""
            try:
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
            self._docks[object_name] = dock
            return dock

        def _build_status_bar(self) -> None:
            status = QStatusBar(self)
            status.setObjectName("workbench.status_bar")
            initial_view = session_state_view(self.session.state.value)
            status.showMessage(
                f"State: {initial_view.marker} {initial_view.label} | Iteration: {self.session.iteration} | Time: {self.session.time:g}"
            )
            self.setStatusBar(status)

        @staticmethod
        def _status_color(status: WorkflowStatus) -> str:
            return {
                WorkflowStatus.COMPLETE: TOKENS["status.success"],
                WorkflowStatus.WARNING: TOKENS["status.warning"],
                WorkflowStatus.ERROR: TOKENS["status.error"],
                WorkflowStatus.NOT_CONFIGURED: TOKENS["text.muted"],
                WorkflowStatus.CONDITIONAL: TOKENS["status.info"],
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
            view = session_state_view(state.simulation_state.value)
            self.statusBar().showMessage(
                f"State: {view.marker} {view.label} | Iteration: {state.execution.iteration} | Time: {state.execution.time:g}"
            )
            self._refresh_properties()
            if hasattr(self, "setup_panel"):
                self.setup_panel.set_diagnostics(state.diagnostics)
            self._refresh_action_availability()

        def _selection_changed(self) -> None:
            items = self.workflow_tree.selectedItems()
            if not items:
                return
            item = items[0]
            stable_id = item.data(0, Qt.ItemDataRole.UserRole)
            kind = item.data(0, Qt.ItemDataRole.UserRole + 1)
            label = item.text(0).rsplit("  [", 1)[0]
            self.application.select(stable_id, kind, label)

        def _selection_state_changed(self, event) -> None:
            """Keep the typed inspector in sync with the selection context."""
            self._refresh_properties()

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
            self._settings.setValue("workspace", self._workspace)
            super().closeEvent(event)

else:

    class CFDXWorkbenchWindow:
        """Clear failure when optional GUI dependencies are unavailable."""

        def __init__(self, *args, **kwargs) -> None:
            raise RuntimeError("PySide6 is required for the CFDX GUI")
