"""Qt results browser backed by Application.results state."""
from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QComboBox, QLabel, QListWidget, QPushButton, QVBoxLayout, QWidget

from .application import Application, ResultsChanged


class ResultsPanel(QWidget):
    """Browse result frames and fields without owning result data."""

    def __init__(self, application: Application, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.application = application
        self.on_open: Callable[[], None] | None = None
        self.setObjectName("workbench.results_panel")
        self.open_button = QPushButton("Open Results Directory…")
        self.frame_list = QListWidget()
        self.field = QComboBox()
        self.status = QLabel("No results loaded")
        layout = QVBoxLayout(self)
        layout.addWidget(self.open_button)
        layout.addWidget(QLabel("Frames"))
        layout.addWidget(self.frame_list, 1)
        layout.addWidget(QLabel("Field"))
        layout.addWidget(self.field)
        layout.addWidget(self.status)
        self.open_button.clicked.connect(self._open_requested)
        self.frame_list.currentItemChanged.connect(self._frame_changed)
        self.field.currentTextChanged.connect(self._field_changed)
        application.events.subscribe(ResultsChanged, self._results_changed)
        self._refresh(application.state)

    def _open_requested(self) -> None:
        if self.on_open:
            self.on_open()

    def _results_changed(self, event) -> None:
        self._refresh(event.state)

    def _refresh(self, state) -> None:
        results = state.results
        self.frame_list.blockSignals(True)
        self.frame_list.clear()
        selected_row = -1
        for row, frame in enumerate(results.frames):
            time = f"t={frame.time:g}" if frame.time is not None else "time=?"
            suffix = " [incomplete]" if not frame.complete else ""
            self.frame_list.addItem(f"{row + 1}: {time}{suffix} — {frame.path}")
            list_item = self.frame_list.item(row)
            list_item.setData(Qt.ItemDataRole.UserRole, frame.stable_id)
            if frame.stable_id == results.selected_frame_id:
                selected_row = row
        if selected_row >= 0:
            self.frame_list.setCurrentRow(selected_row)
        self.frame_list.blockSignals(False)
        self.field.blockSignals(True)
        self.field.clear()
        self.field.addItems(results.field_names)
        if results.selected_field:
            self.field.setCurrentText(results.selected_field)
        self.field.blockSignals(False)
        self.status.setText(
            f"{len(results.frames)} frame(s)"
            + (f" | {results.selected_frame_id} | {results.selected_field}" if results.selected_frame_id else "")
        )

    def _frame_changed(self, current, _previous) -> None:
        if current is not None:
            self.application.select_result_frame(current.data(Qt.ItemDataRole.UserRole))

    def _field_changed(self, field: str) -> None:
        if field:
            self.application.select_result_field(field)
