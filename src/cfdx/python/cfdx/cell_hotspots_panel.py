"""Optional Qt viewer for native accepted-iteration per-cell diagnostics.

The sparse top-K overlay is never presented as a full-cell heatmap or a
physical convergence certificate. A bounded trace is parsed on a worker thread.
"""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QObject, QRunnable, QThreadPool, QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .cell_convergence import CellConvergenceError, read_cell_trace
from .cell_field_map import CellFieldMap, CellMapError


class _ReadSignals(QObject):
    finished = Signal(int, object, object)


class _ReadJob(QRunnable):
    def __init__(self, revision: int, path: Path) -> None:
        super().__init__()
        self.revision = revision
        self.path = path
        self.events = _ReadSignals()

    def run(self) -> None:
        try:
            trace = read_cell_trace(self.path)
            error = None
        except (CellConvergenceError, OSError, ValueError) as exc:
            trace = None
            error = str(exc)
        self.events.finished.emit(self.revision, trace, error)


class _MapJob(QRunnable):
    def __init__(self, revision, source, iteration, metric, mesh, n_cells):
        super().__init__()
        self.revision, self.source = revision, source
        self.iteration, self.metric = iteration, metric
        self.mesh, self.n_cells = mesh, n_cells
        self.events = _ReadSignals()

    def run(self):
        try:
            diagnostic = CellFieldMap(self.source)
            if (
                Path(diagnostic.mesh_source).resolve() != Path(self.mesh).resolve()
                or diagnostic.n_cells != self.n_cells
            ):
                raise CellMapError(
                    "full map mesh is incompatible with selected top-K trace"
                )
            values, clim, title = diagnostic.scalar(self.metric, self.iteration)
            payload = {
                "mesh_source": diagnostic.mesh_source,
                "n_cells": diagnostic.n_cells,
                "iteration": self.iteration,
                "values": values,
                "clim": clim,
                "title": title,
            }
            diagnostic.close()
            self.events.finished.emit(self.revision, payload, None)
        except (CellMapError, OSError, ValueError, RuntimeError) as exc:
            self.events.finished.emit(self.revision, None, str(exc))


class CellHotspotsPanel(QWidget):
    """Progressive spatial diagnosis: inspect, compare and locate actual cell IDs."""

    frame_selected = Signal(object, object, str)  # trace, frame, metric
    full_map_selected = Signal(object)  # values + scale, on originating mesh
    cell_selected = Signal(object, object)  # trace, sample

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("workbench.cell_hotspots")
        self._source: Path | None = None
        self._signature: tuple[int, int] | None = None
        self._revision = 0
        self._busy = False
        self._map_revision = 0
        self._map_jobs = {}
        self._read_jobs: dict[int, _ReadJob] = {}
        self.trace = None
        self.status = QLabel(
            "Enable spatial diagnostics before Run, or open a real cell_convergence.jsonl. "
            "Only sampled hotspots are shown; no full-field convergence claim."
        )
        self.status.setWordWrap(True)
        self.status.setObjectName("workbench.cell_hotspots.status")
        self.open_button = QPushButton("Open trace…", self)
        self.open_button.setObjectName("workbench.cell_hotspots.open")
        self.open_button.clicked.connect(self._choose_file)
        self.reload_button = QPushButton("Refresh", self)
        self.reload_button.setObjectName("workbench.cell_hotspots.refresh")
        self.reload_button.clicked.connect(lambda: self.refresh(force=True))
        self.export_button = QPushButton("Export report", self)
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self._export)
        buttons = QHBoxLayout()
        buttons.addWidget(self.open_button)
        buttons.addWidget(self.reload_button)
        buttons.addWidget(self.export_button)
        self.frame_choice = QComboBox(self)
        self.frame_choice.setObjectName("workbench.cell_hotspots.iteration")
        self.frame_choice.currentIndexChanged.connect(self._show_frame)
        self.metric_choice = QComboBox(self)
        self.metric_choice.setObjectName("workbench.cell_hotspots.metric")
        self.metric_choice.addItem("Momentum row defect / RHS", "momentum")
        self.metric_choice.addItem("Mass imbalance (absolute)", "mass")
        self.metric_choice.addItem("Non-decrease streak", "persistence")
        self.metric_choice.currentIndexChanged.connect(self._show_frame)
        choices = QHBoxLayout()
        choices.addWidget(self.frame_choice)
        choices.addWidget(self.metric_choice)
        self.table = QTableWidget(0, 7, self)
        self.table.setObjectName("workbench.cell_hotspots.table")
        self.table.setHorizontalHeaderLabels(
            ("Cell", "Momentum / RHS", "|mass defect|", "Non-decrease", "X", "Y", "Z")
        )
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._select_cell)
        self.history_note = QLabel(
            "Select a cell to see only its recorded top-K snapshots."
        )
        self.history_note.setWordWrap(True)
        self.history_table = QTableWidget(0, 4, self)
        self.history_table.setObjectName("workbench.cell_hotspots.history")
        self.history_table.setHorizontalHeaderLabels(
            ("Accepted iter", "Momentum / RHS", "|mass defect|", "Non-decrease")
        )
        self.history_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.history_table.horizontalHeader().setStretchLastSection(True)
        self.history_table.setMaximumHeight(150)
        self.notice = QLabel(
            "Residual is A(U,p)U − b(U,p) reassembled at the accepted state. "
            "Momentum values use per-component global RHS scales. "
            "Three non-decreasing transitions are a diagnostic flag, not proof of stagnation. "
            "When a full cell-capable dataset is loaded, sampled cells are coloured on the mesh "
            "with a blue→red scale; unsampled cells remain neutral grey."
        )
        self.notice.setWordWrap(True)
        self.notice.setObjectName("workbench.cell_hotspots.disclaimer")
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addLayout(buttons)
        layout.addLayout(choices)
        layout.addWidget(self.table)
        layout.addWidget(self.history_note)
        layout.addWidget(self.history_table)
        layout.addWidget(self.notice)
        self._pool = QThreadPool.globalInstance()
        self._watch = QTimer(self)
        self._watch.setInterval(2500)
        self._watch.timeout.connect(self.refresh)
        self._watch.start()

    @property
    def source(self) -> Path | None:
        return self._source

    def set_source(self, path: str | Path | None, *, automatic: bool = False) -> None:
        source = Path(path) if path is not None else None
        if source == self._source:
            return
        self._revision += 1
        self._map_revision += 1
        self._source = source
        self._signature = None
        self.trace = None
        self.frame_choice.clear()
        self.table.setRowCount(0)
        self.history_table.setRowCount(0)
        self.export_button.setEnabled(False)
        if source is None:
            self.status.setText("No per-cell diagnostic is attached to this run.")
        else:
            self.status.setText(
                f"{'Observing native run' if automatic else 'Reading spatial diagnostic'}: "
                f"{source} (waiting for the first committed record)"
            )
            self.refresh(force=True)

    def _choose_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open CFDX per-cell trace", "", "CFDX JSONL (*.jsonl);;All files (*)"
        )
        if path:
            self.set_source(path)

    def refresh(self, *, force: bool = False) -> None:
        source = self._source
        if source is None or self._busy:
            return
        try:
            info = source.stat()
        except FileNotFoundError:
            # A previous trace at this path is not evidence for a new run.
            if self.trace is not None:
                self.trace = None
                self._signature = None
                self.frame_choice.clear()
                self.table.setRowCount(0)
                self.history_table.setRowCount(0)
                self.export_button.setEnabled(False)
                self.status.setText(
                    "Cell trace no longer exists; waiting for a new native write."
                )
            return  # Solver has not yet written its opt-in header.
        except OSError as exc:
            self.status.setText(f"Cannot inspect trace: {exc}")
            return
        signature = (info.st_size, info.st_mtime_ns)
        if not force and signature == self._signature:
            return
        if info.st_size > 32 * 1024 * 1024:
            self.trace = None
            self.frame_choice.clear()
            self.table.setRowCount(0)
            self.history_table.setRowCount(0)
            self.export_button.setEnabled(False)
            self.status.setText(
                "Trace exceeds 32 MiB: increase --cell-hotspots-every or inspect via CLI."
            )
            return
        self._revision += 1
        revision = self._revision
        self._busy = True
        job = _ReadJob(revision, source)
        self._read_jobs[revision] = job
        job.events.finished.connect(self._read_finished)
        self._pool.start(job)

    def _read_finished(self, revision: int, trace, error: str | None) -> None:
        self._read_jobs.pop(revision, None)
        self._busy = False
        if revision != self._revision or self._source is None:
            return  # The operator selected another file in the meantime.
        if error:
            self.status.setText(
                f"Diagnostic rejected (no unverified data displayed): {error}"
            )
            self.trace = None
            self.frame_choice.clear()
            self.table.setRowCount(0)
            self.history_table.setRowCount(0)
            self.export_button.setEnabled(False)
            return
        self.trace = trace
        try:
            info = self._source.stat()
            self._signature = (info.st_size, info.st_mtime_ns)
        except OSError:
            self._signature = None
        self.export_button.setEnabled(True)
        self.status.setText(
            f"{trace.n_cells:,} mesh cells | {len(trace.frames)} snapshots | "
            f"{'complete stream' if trace.complete else 'live / uncompleted stream'} | "
            f"source: {trace.mesh_source}"
        )
        selected = self.frame_choice.currentData()
        self.frame_choice.blockSignals(True)
        self.frame_choice.clear()
        for frame in trace.frames:
            self.frame_choice.addItem(
                f"Accepted iteration {frame.iteration}", frame.iteration
            )
        if selected is not None:
            found = self.frame_choice.findData(selected)
            if found >= 0:
                self.frame_choice.setCurrentIndex(found)
            else:
                self.frame_choice.setCurrentIndex(self.frame_choice.count() - 1)
        else:
            self.frame_choice.setCurrentIndex(self.frame_choice.count() - 1)
        self.frame_choice.blockSignals(False)
        self._show_frame()

    def _current_frame(self):
        if self.trace is None:
            return None
        iteration = self.frame_choice.currentData()
        return next(
            (frame for frame in self.trace.frames if frame.iteration == iteration), None
        )

    def _show_frame(self, *_args) -> None:
        frame = self._current_frame()
        trace = self.trace
        self.table.setRowCount(0)
        self.history_table.setRowCount(0)
        self.history_note.setText(
            "Select a cell to see its actually recorded samples (gaps are unknown)."
        )
        if frame is None or trace is None:
            return
        metric = self.metric_choice.currentData() or "momentum"
        ranked = frame.ranked(metric, 72)
        self.table.setRowCount(len(ranked))
        for row, cell in enumerate(ranked):
            columns = (
                cell.cell_id,
                f"{cell.momentum_scaled:.5g}",
                f"{cell.mass_imbalance_abs:.5g}",
                cell.nondecrease_streak,
                *(f"{coord:.6g}" for coord in cell.center),
            )
            for col, value in enumerate(columns):
                item = QTableWidgetItem(str(value))
                if col == 0:
                    item.setData(0x0100, cell)  # Qt.UserRole
                self.table.setItem(row, col, item)
        self.frame_selected.emit(trace, frame, metric)
        self._map_revision += 1
        if self._source is not None:
            native_map = self._source.with_name("cell_convergence_map.bin")
            if native_map.is_file():
                revision = self._map_revision
                job = _MapJob(
                    revision,
                    native_map,
                    frame.iteration,
                    metric,
                    trace.mesh_source,
                    trace.n_cells,
                )
                self._map_jobs[revision] = job
                job.events.finished.connect(self._map_finished)
                self._pool.start(job)

    def _map_finished(self, revision, payload, error) -> None:
        self._map_jobs.pop(revision, None)
        if revision != self._map_revision or self.trace is None:
            return
        if error:
            self.status.setText(
                "Full-cell map unavailable; sparse hotspots remain: " + error
            )
        else:
            self.full_map_selected.emit(payload)

    def _select_cell(self) -> None:
        row = self.table.currentRow()
        if row < 0 or self.trace is None:
            return
        item = self.table.item(row, 0)
        cell = item.data(0x0100) if item is not None else None
        if cell is not None:
            history = self.trace.observed_history(cell.cell_id)
            self.history_table.setRowCount(len(history))
            for line, (accepted_iter, sample) in enumerate(history):
                values = (
                    accepted_iter,
                    f"{sample.momentum_scaled:.5g}",
                    f"{sample.mass_imbalance_abs:.5g}",
                    sample.nondecrease_streak,
                )
                for column, value in enumerate(values):
                    self.history_table.setItem(
                        line, column, QTableWidgetItem(str(value))
                    )
            self.history_note.setText(
                f"Cell {cell.cell_id} • {len(history)} recorded snapshots. "
                "Missing snapshots are unknown; the streak itself was computed by the solver "
                "over every accepted nonlinear iteration."
            )
            self.cell_selected.emit(self.trace, cell)

    def _export(self) -> None:
        if self.trace is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export CFDX cell hotspot diagnostic",
            "cell_hotspots_report.json",
            "JSON (*.json)",
        )
        if not path:
            return
        try:
            Path(path).write_text(
                json.dumps(
                    self.trace.to_dict(), indent=2, sort_keys=True, allow_nan=False
                )
                + "\n",
                encoding="utf-8",
            )
        except (OSError, ValueError) as exc:
            self.status.setText(f"Cannot export diagnostic report: {exc}")
