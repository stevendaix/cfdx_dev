"""Qt adapter for the headless Application Run Center."""
from __future__ import annotations

from .application import Application, ApplicationStateChanged

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QTextEdit, QVBoxLayout, QWidget


class RunCenterPanel(QWidget):
    """Controls and live output for the application execution controller."""

    def __init__(self, application: Application, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.application = application
        self.setObjectName("workbench.run_center")
        self.status = QLabel("No execution attached")
        self.metrics = QLabel("Iteration: — | Time: — | CFL: —")
        self.output = QTextEdit()
        self.output.setReadOnly(True)
        self.output.setPlaceholderText("Solver output will appear here")
        self.run_button = QPushButton("Run")
        self.pause_button = QPushButton("Pause")
        self.resume_button = QPushButton("Resume")
        self.stop_button = QPushButton("Stop")
        controls = QHBoxLayout()
        for button in (self.run_button, self.pause_button, self.resume_button, self.stop_button):
            controls.addWidget(button)
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.metrics)
        layout.addLayout(controls)
        layout.addWidget(self.output)
        self.run_button.clicked.connect(self._run)
        self.pause_button.clicked.connect(self._pause)
        self.resume_button.clicked.connect(self._resume)
        self.stop_button.clicked.connect(self._stop)
        application.events.subscribe(ApplicationStateChanged, self._state_changed)
        self._refresh(application.state)

    def _run(self) -> None:
        self.application.run()

    def _pause(self) -> None:
        self.application.pause()

    def _resume(self) -> None:
        self.application.resume()

    def _stop(self) -> None:
        self.application.stop()

    def _state_changed(self, event) -> None:
        self._refresh(event.state)

    def _refresh(self, state) -> None:
        center = self.application.run_center
        if center is None:
            self.status.setText("No execution controller attached")
            self.metrics.setText("Iteration: — | Time: — | CFL: —")
            self.output.clear()
            for button in (self.run_button, self.pause_button, self.resume_button, self.stop_button):
                button.setEnabled(False)
            return
        center_state = center.state
        self.status.setText(f"State: {center_state.simulation_state.value}")
        cfl = center_state.latest_metrics.cfl if center_state.latest_metrics else None
        cfl_text = f"{cfl:g}" if cfl is not None else "—"
        self.metrics.setText(
            f"Iteration: {center_state.iteration} | Time: {center_state.time:g} | CFL: {cfl_text}"
        )
        self.output.setPlainText("\n".join(center_state.output))
        self.output.verticalScrollBar().setValue(self.output.verticalScrollBar().maximum())
        running = center_state.simulation_state.value == "RUNNING"
        paused = center_state.simulation_state.value == "PAUSED"
        self.run_button.setEnabled(not running and not paused)
        self.pause_button.setEnabled(running)
        self.resume_button.setEnabled(paused)
        self.stop_button.setEnabled(running or paused)
