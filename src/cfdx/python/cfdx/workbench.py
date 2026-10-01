"""Dockable CFDX Workbench shell.

The shell deliberately contains presentation placeholders only. Domain actions
remain owned by the application/session layer and the existing GUI adapters.
"""
from __future__ import annotations

from .session import CFDXSession

try:
    from PySide6.QtCore import QSettings, Qt
    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import (
        QDockWidget,
        QLabel,
        QMainWindow,
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

    class CFDXWorkbenchWindow(QMainWindow):
        """Initial Workbench composition root with stable dock object names."""

        SETTINGS_ORGANIZATION = "CFDX"
        SETTINGS_APPLICATION = "Workbench"

        def __init__(self, session: CFDXSession | None = None) -> None:
            super().__init__()
            self.session = session or CFDXSession()
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
            setup = QTreeWidgetItem(["SETUP"])
            for label in ("Geometry", "Mesh", "Physics", "Materials", "Boundaries", "Numerics", "Solver"):
                setup.addChild(QTreeWidgetItem([label]))
            run = QTreeWidgetItem(["RUN"])
            for label in ("Check", "Initialize", "Run", "Monitor", "Checkpoint"):
                run.addChild(QTreeWidgetItem([label]))
            results = QTreeWidgetItem(["RESULTS"])
            for label in ("Fields", "Contours", "Vectors", "Slices", "Probes", "Reports"):
                results.addChild(QTreeWidgetItem([label]))
            self.workflow_tree.addTopLevelItems([setup, run, results])
            setup.setExpanded(True)
            run.setExpanded(True)
            results.setExpanded(True)
            self._add_dock("Workflow", "workbench.dock.workflow", self.workflow_tree, Qt.DockWidgetArea.LeftDockWidgetArea)

            viewport = QLabel("3D VIEWPORT\n\nRenderer adapter placeholder")
            viewport.setObjectName("workbench.viewport")
            viewport.setAlignment(Qt.AlignmentFlag.AlignCenter)
            viewport.setMinimumSize(480, 320)
            self._add_dock("Viewport", "workbench.dock.viewport", viewport, Qt.DockWidgetArea.RightDockWidgetArea)
            self.setCentralWidget(viewport)

            properties = QTextEdit()
            properties.setObjectName("workbench.properties")
            properties.setReadOnly(True)
            properties.setPlaceholderText("Select an item to inspect its properties")
            self._add_dock("Properties", "workbench.dock.properties", properties, Qt.DockWidgetArea.RightDockWidgetArea)

            monitor = QTextEdit()
            monitor.setObjectName("workbench.monitor")
            monitor.setReadOnly(True)
            monitor.setPlainText("MONITORS\nNo execution attached")
            self._add_dock("Monitors / Console", "workbench.dock.monitor", monitor, Qt.DockWidgetArea.BottomDockWidgetArea)

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
