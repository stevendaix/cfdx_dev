"""Dockable CFDX Workbench shell.

The shell deliberately contains presentation placeholders only. Domain actions
remain owned by the application/session layer and the existing GUI adapters.
"""
from __future__ import annotations

from .application import Application, ApplicationStateChanged, WorkflowStatus, workflow_children
from .session import CFDXSession

try:
    from PySide6.QtCore import QSettings, Qt
    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import (
        QDockWidget,
        QFileDialog,
        QFormLayout,
        QLabel,
        QLineEdit,
        QMainWindow,
        QPushButton,
        QStatusBar,
        QTextEdit,
        QTreeWidget,
        QTreeWidgetItem,
        QToolBar,
        QWidget,
    )
except ImportError:  # pragma: no cover
    QMainWindow = object


if QMainWindow is not object:
    from .run_center_panel import RunCenterPanel
    from .results_panel import ResultsPanel

    class CFDXWorkbenchWindow(QMainWindow):
        """Initial Workbench composition root with stable dock object names."""

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
            self.application.events.subscribe(ApplicationStateChanged, self._state_changed)
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

            for label in ("New Project", "Open", "Save"):
                action = QAction(label, self)
                action.setObjectName(f"workbench.action.{label.lower().replace(' ', '_')}")
                action.setEnabled(False)
                toolbar.addAction(action)
            toolbar.addSeparator()
            for label in ("Check", "Run", "Stop"):
                action = QAction(label, self)
                action.setObjectName(f"workbench.action.{label.lower()}")
                action.setEnabled(False)
                toolbar.addAction(action)
            toolbar.addSeparator()
            search = QAction("Search", self)
            search.setObjectName("workbench.action.search")
            search.setEnabled(False)
            toolbar.addAction(search)
            self.toolbar = toolbar

        def _build_docks(self) -> None:
            self.workflow_tree = QTreeWidget()
            self.workflow_tree.setObjectName("workbench.workflow_tree")
            self.workflow_tree.setHeaderLabel("Workflow")
            self.workflow_tree.itemSelectionChanged.connect(self._selection_changed)
            self._add_dock("Workflow", "workbench.dock.workflow", self.workflow_tree, Qt.DockWidgetArea.LeftDockWidgetArea)
            self._refresh_workflow(self._application_state)
            self.results_panel = ResultsPanel(self.application)
            self.results_panel.on_open = self._open_results_directory
            self._add_dock("Results", "workbench.dock.results", self.results_panel, Qt.DockWidgetArea.LeftDockWidgetArea)

            viewport = QLabel("3D VIEWPORT\n\nRenderer adapter placeholder")
            viewport.setObjectName("workbench.viewport")
            viewport.setAlignment(Qt.AlignmentFlag.AlignCenter)
            viewport.setMinimumSize(480, 320)
            self.setCentralWidget(viewport)

            properties = QTextEdit()
            properties = self._build_properties_panel()
            self._add_dock("Properties", "workbench.dock.properties", properties, Qt.DockWidgetArea.RightDockWidgetArea)

            monitor = RunCenterPanel(self.application)
            monitor.setObjectName("workbench.monitor")
            self._add_dock("Monitors / Console", "workbench.dock.monitor", monitor, Qt.DockWidgetArea.BottomDockWidgetArea)

        def _build_properties_panel(self) -> QWidget:
            panel = QWidget()
            panel.setObjectName("workbench.properties")
            panel._form = QFormLayout(panel)
            panel._form.addRow(QLabel("Select a workflow item"))
            self.properties_panel = panel
            self._refresh_properties()
            return panel

        def _open_results_directory(self) -> None:
            directory = QFileDialog.getExistingDirectory(self, "Open CFDX Results Directory")
            if not directory:
                return
            try:
                self.application.open_results(directory)
            except (OSError, ValueError) as exc:
                self.statusBar().showMessage(f"Results error: {exc}")

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
                button.clicked.connect(lambda _checked=False, p=prop, e=editor: self._apply_property(p.key, e.text()))
                form.addRow(prop.label, editor)
                form.addRow("", button)

        def _apply_property(self, key: str, raw_value: str) -> None:
            try:
                value = float(raw_value) if key == "numerics.cfl" else raw_value
                self.application.set_property(key, value)
            except (KeyError, TypeError, ValueError) as exc:
                self.statusBar().showMessage(f"Property error: {exc}")
            else:
                self._refresh_properties()

        def _add_dock(self, title: str, object_name: str, widget: QWidget, area: Qt.DockWidgetArea) -> QDockWidget:
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
        def _status_marker(status: WorkflowStatus) -> str:
            return {
                WorkflowStatus.COMPLETE: "[OK]",
                WorkflowStatus.WARNING: "[!]",
                WorkflowStatus.ERROR: "[X]",
                WorkflowStatus.NOT_CONFIGURED: "[ ]",
            }[status]

        def _refresh_workflow(self, state) -> None:
            self.workflow_tree.clear()
            roots = workflow_children(state.workflow, None)
            for step in roots:
                root = QTreeWidgetItem([f"{step.label}  {self._status_marker(step.status)}"])
                root.setData(0, Qt.ItemDataRole.UserRole, step.id)
                root.setData(0, Qt.ItemDataRole.UserRole + 1, step.kind)
                root.setToolTip(0, step.message or step.status.value)
                self.workflow_tree.addTopLevelItem(root)
                for child in workflow_children(state.workflow, step.id):
                    item = QTreeWidgetItem([f"{child.label}  {self._status_marker(child.status)}"])
                    item.setData(0, Qt.ItemDataRole.UserRole, child.id)
                    item.setData(0, Qt.ItemDataRole.UserRole + 1, child.kind)
                    item.setToolTip(0, child.message or child.status.value)
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
