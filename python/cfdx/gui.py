"""Optional PySide6 GUI shell backed by the headless CFDX session."""
from __future__ import annotations

from .session import CFDXSession, SimulationState
from .tui import TuiRenderer

try:
    from PySide6.QtCore import QObject, Qt, Signal
    from PySide6.QtWidgets import (
        QApplication, QHBoxLayout, QLabel, QMainWindow, QPlainTextEdit,
        QPushButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
    )
except ImportError:  # pragma: no cover
    QApplication = None
    QObject = object
    Signal = None


if QApplication is not None:
    class SessionSignals(QObject):
        state_changed = Signal(object)
        output = Signal(str, bool)

    class CFDXMainWindow(QMainWindow):
        """Minimal Fluent-like shell with shared session controls."""

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
            self.run_button = QPushButton("Run")
            self.pause_button = QPushButton("Pause")
            self.stop_button = QPushButton("Stop")
            controls = QHBoxLayout()
            controls.addWidget(self.run_button)
            controls.addWidget(self.pause_button)
            controls.addWidget(self.stop_button)
            right = QVBoxLayout()
            right.addWidget(self.status)
            right.addLayout(controls)
            right.addWidget(self.log)
            layout.addWidget(self.tree, 1)
            container = QWidget()
            container.setLayout(right)
            layout.addWidget(container, 3)
            self.setCentralWidget(central)

            self.signals.state_changed.connect(self._refresh_status)
            self.signals.output.connect(self.append_output)
            self.run_button.clicked.connect(self._run)
            self.pause_button.clicked.connect(self._pause)
            self.stop_button.clicked.connect(self._stop)
            self.refresh()

        def refresh(self) -> None:
            self.tree.clear()
            for node in self.session.case_tree():
                item = QTreeWidgetItem([node.label])
                item.setData(0, Qt.ItemDataRole.UserRole, node.id)
                self.tree.addTopLevelItem(item)
            self._refresh_status(self.session.state)

        def _refresh_status(self, state: SimulationState) -> None:
            self.status.setText(
                f"State: {state.value} | Iteration: {self.session.iteration} | "
                f"Time: {self.session.time:g}"
            )
            self.run_button.setEnabled(state in {
                SimulationState.CREATED, SimulationState.READY,
                SimulationState.PAUSED, SimulationState.STOPPED,
            })
            self.pause_button.setEnabled(state is SimulationState.RUNNING)
            self.stop_button.setEnabled(state in {
                SimulationState.RUNNING, SimulationState.PAUSED,
            })

        def _run(self) -> None:
            self.session.run()
            self.signals.state_changed.emit(self.session.state)

        def _pause(self) -> None:
            self.session.pause()
            self.signals.state_changed.emit(self.session.state)

        def _stop(self) -> None:
            self.session.stop()
            self.signals.state_changed.emit(self.session.state)

        def append_output(self, line: str, is_stderr: bool = False) -> None:
            self.log.appendPlainText(("[stderr] " if is_stderr else "") + line)

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
        def __init__(self, session: CFDXSession | None = None) -> None:
            raise RuntimeError("PySide6 is required for the CFDX GUI")

    def create_application(argv: list[str] | None = None):
        raise RuntimeError("PySide6 is required for the CFDX GUI")

    def launch(session: CFDXSession | None = None, argv: list[str] | None = None) -> int:
        raise RuntimeError("PySide6 is required for the CFDX GUI")
