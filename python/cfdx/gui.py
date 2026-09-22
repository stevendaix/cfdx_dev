"""Optional PySide6 GUI shell backed by the headless CFDX session."""
from __future__ import annotations

from typing import Callable

from .session import CFDXSession, SimulationState
from .tui import TuiRenderer

try:
    from PySide6.QtCore import QObject, Signal
    from PySide6.QtWidgets import (
        QApplication,
        QHBoxLayout,
        QLabel,
        QMainWindow,
        QPlainTextEdit,
        QTreeWidget,
        QTreeWidgetItem,
        QVBoxLayout,
        QWidget,
    )
except ImportError:  # pragma: no cover - exercised through optional-dependency tests
    QApplication = None
    QObject = object
    Signal = None


if QApplication is not None:
    class SessionSignals(QObject):
        """Qt bridge for state and log events."""

        state_changed = Signal(object)
        output = Signal(str, bool)

    class CFDXMainWindow(QMainWindow):
        """Minimal Fluent-like shell: case tree, status and TUI/log pane."""

        def __init__(self, session: CFDXSession | None = None) -> None:
            super().__init__()
            self.session = session or CFDXSession()
            self.signals = SessionSignals()
            self.setWindowTitle(f"CFDX — {self.session.case.name}")
            self.resize(1000, 650)

            central = QWidget()
            layout = QHBoxLayout(central)
            self.tree = QTreeWidget()
            self.tree.setHeaderLabels(["Case"])
            self.log = QPlainTextEdit()
            self.log.setReadOnly(True)
            self.status = QLabel()
            right = QVBoxLayout()
            right.addWidget(self.status)
            right.addWidget(self.log)
            layout.addWidget(self.tree, 1)
            container = QWidget()
            container.setLayout(right)
            layout.addWidget(container, 3)
            self.setCentralWidget(central)

            self.signals.state_changed.connect(self._refresh_status)
            self.signals.output.connect(self.append_output)
            self.refresh()

        def refresh(self) -> None:
            self.tree.clear()
            for node in self.session.case_tree():
                item = QTreeWidgetItem([node.label])
                item.setData(0, 32, node.id)
                self.tree.addTopLevelItem(item)
            self._refresh_status(self.session.state)

        def _refresh_status(self, state: SimulationState) -> None:
            self.status.setText(
                f"State: {state.value} | Iteration: {self.session.iteration} | "
                f"Time: {self.session.time:g}"
            )

        def append_output(self, line: str, is_stderr: bool = False) -> None:
            prefix = "[stderr] " if is_stderr else ""
            self.log.appendPlainText(prefix + line)

        def render_tui(self) -> str:
            return TuiRenderer.render(self.session)


    def create_application(argv: list[str] | None = None) -> QApplication:
        return QApplication.instance() or QApplication(argv or [])


    def launch(session: CFDXSession | None = None, argv: list[str] | None = None) -> int:
        app = create_application(argv)
        window = CFDXMainWindow(session)
        window.show()
        return app.exec()
else:
    class CFDXMainWindow:
        """Placeholder that gives a clear error when GUI dependencies are absent."""

        def __init__(self, session: CFDXSession | None = None) -> None:
            raise RuntimeError("PySide6 is required for the CFDX GUI")


    def create_application(argv: list[str] | None = None):
        raise RuntimeError("PySide6 is required for the CFDX GUI")


    def launch(session: CFDXSession | None = None, argv: list[str] | None = None) -> int:
        raise RuntimeError("PySide6 is required for the CFDX GUI")
