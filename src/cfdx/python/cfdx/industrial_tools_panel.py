"""Industrial Workbench tools: bounded data inspection and verified comparisons."""

from __future__ import annotations

import json

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal
from PySide6.QtWidgets import (
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .conservative_remap import conservative_box_remap
from .convex_polyhedral_remap import conservative_convex_remap
from .progressive_fields import stream_field_statistics
from .simulation_compare import compare_run_reports, compare_scalar_fields


class _Events(QObject):
    done = Signal(object, object)


class _Job(QRunnable):
    def __init__(self, callback):
        super().__init__()
        self.callback = callback
        self.events = _Events()

    def run(self):
        try:
            self.events.done.emit(self.callback(), None)
        except Exception as exc:
            self.events.done.emit(None, str(exc))


class IndustrialToolsPanel(QWidget):
    """Docks into Analyze; scientific computations never run in the GUI thread."""

    show_difference = Signal(str, str, str)
    volume_partition_requested = Signal(str, int, int, object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("workbench.industrial_tools")
        self._pool = QThreadPool.globalInstance()
        self._jobs = []
        self._last_report = None
        layout = QVBoxLayout(self)
        self.tabs = QTabWidget(self)
        layout.addWidget(self.tabs)
        progressive = QWidget()
        form = QFormLayout(progressive)
        self.source_file = QLineEdit()
        self.field = QLineEdit("/fields/scalar/p")
        self.stats_text = QPlainTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setPlaceholderText("Load an HDF5 field using bounded chunks.")
        row = QHBoxLayout()
        select = QPushButton("Open HDF5…")
        select.clicked.connect(lambda: self._choose(self.source_file))
        inspect = QPushButton("Inspect progressive field")
        inspect.clicked.connect(self.inspect)
        row.addWidget(select)
        row.addWidget(inspect)
        form.addRow("File", self.source_file)
        form.addRow("Scalar dataset", self.field)
        form.addRow(row)
        form.addRow(self.stats_text)
        self.tabs.addTab(progressive, "Large-field inspector")

        tile_tab = QWidget(self)
        vf = QFormLayout(tile_tab)
        self.volume_file = QLineEdit()
        self.volume_file.setPlaceholderText("Source native CFDX mesh HDF5")
        self.volume_start = QSpinBox()
        self.volume_start.setRange(0, 1_000_000_000)
        self.volume_count = QSpinBox()
        self.volume_count.setRange(1, 512)
        self.volume_count.setValue(128)
        self.volume_scalar = QLineEdit()
        self.volume_scalar.setPlaceholderText("Optional /fields/scalar/T")
        vf.addRow("Volume HDF5", self.volume_file)
        vf.addRow("First original cell ID", self.volume_start)
        vf.addRow("Cells per tile (1..512)", self.volume_count)
        vf.addRow("Optional cell scalar", self.volume_scalar)
        tile_button = QPushButton("Load one actual polyhedral volume tile")
        tile_button.setToolTip(
            "I/O off Qt thread; exact cell IDs; VTK CPU/OpenGL renderer; not GPU compute"
        )
        tile_button.clicked.connect(self.request_volume_partition)
        vf.addRow(tile_button)
        self.volume_note = QLabel(
            "Only the selected cell range is loaded. Other cells are absent, not converged or empty. "
        )
        self.volume_note.setWordWrap(True)
        vf.addRow(self.volume_note)
        self.tabs.addTab(tile_tab, "Volume tiles (CPU/VTK)")

        comparison = QWidget()
        cf = QFormLayout(comparison)
        self.left_file, self.right_file = QLineEdit(), QLineEdit()
        self.compare_field = QLineEdit("/fields/scalar/p")
        self.compare_text = QPlainTextEdit()
        self.compare_text.setReadOnly(True)
        for title, control in (
            ("Reference HDF5", self.left_file),
            ("Candidate HDF5", self.right_file),
        ):
            row = QHBoxLayout()
            row.addWidget(control)
            button = QPushButton("Browse…")
            button.clicked.connect(
                lambda _checked=False, target=control: self._choose(target)
            )
            row.addWidget(button)
            cf.addRow(title, row)
        cf.addRow("Field", self.compare_field)
        compare_button = QPushButton("Compare compatible simulations")
        compare_button.clicked.connect(self.compare)
        cf.addRow(compare_button)
        self.difference_button = QPushButton("Show signed difference in 3D (bounded)")
        self.difference_button.setEnabled(False)
        self.difference_button.clicked.connect(
            lambda: self.show_difference.emit(
                self.left_file.text(), self.right_file.text(), self.compare_field.text()
            )
        )
        cf.addRow(self.difference_button)
        cf.addRow(self.compare_text)
        remap_button = QPushButton("Conservative comparison (Cartesian hexahedra)")
        remap_button.setToolTip(
            "Exact volume overlap; only validated axis-aligned 8-vertex cells. "
            "Refuses curved/nonmatching domains and unbounded meshes."
        )
        remap_button.clicked.connect(self.compare_conservative)
        cf.addRow(remap_button)
        convex_button = QPushButton(
            "Conservative comparison (verified convex polyhedra)"
        )
        convex_button.setToolTip(
            "First-order numerical volume overlap via SciPy Qhull with conservative checks. "
            "Rejects concave/warped cells, incompatible domains and large problems."
        )
        convex_button.clicked.connect(self.compare_convex)
        cf.addRow(convex_button)
        self.tabs.addTab(comparison, "Compare runs")

        monitor_panel = QWidget(self)
        mf = QFormLayout(monitor_panel)
        self.reference_run, self.candidate_run = QLineEdit(), QLineEdit()
        mf.addRow("Reference output directory", self.reference_run)
        mf.addRow("Candidate output directory", self.candidate_run)
        telemetry_button = QPushButton(
            "Compare recorded solver metrics and wall forces"
        )
        telemetry_button.clicked.connect(self.compare_telemetry)
        mf.addRow(telemetry_button)
        self.telemetry_text = QPlainTextEdit()
        self.telemetry_text.setReadOnly(True)
        mf.addRow(self.telemetry_text)
        self.tabs.addTab(monitor_panel, "Solver metrics")
        self.disclaimer = QLabel(
            "Data diagnostics only: identical mesh, unit and time checks are mandatory. Conservative remapping is opt-in and restricted to Cartesian boxes; never a V&V verdict."
        )
        self.disclaimer.setWordWrap(True)
        layout.addWidget(self.disclaimer)

    def _choose(self, target):
        chosen, _ = QFileDialog.getOpenFileName(
            self, "Choose CFDX HDF5", "", "CFDX HDF5 (*.h5);;All files (*)"
        )
        if chosen:
            target.setText(chosen)
            self.difference_button.setEnabled(False)

    def _launch(self, callback, target, on_success=None):
        target.setPlainText("Analyzing... (data worker, GUI remains responsive)")
        job = _Job(callback)
        self._jobs.append(job)

        def finished(result, error):
            if job in self._jobs:
                self._jobs.remove(job)
            if error:
                target.setPlainText("Analysis rejected: " + error)
                self.difference_button.setEnabled(False)
            else:
                if on_success is not None:
                    on_success(result)
                target.setPlainText(
                    json.dumps(result, indent=2, allow_nan=False, ensure_ascii=False)
                )

        job.events.done.connect(finished)
        self._pool.start(job)

    def inspect(self):
        from dataclasses import asdict

        path, field = self.source_file.text(), self.field.text()
        self._launch(
            lambda: asdict(stream_field_statistics(path, field)), self.stats_text
        )

    def compare(self):
        left, right, field = (
            self.left_file.text(),
            self.right_file.text(),
            self.compare_field.text(),
        )
        self.difference_button.setEnabled(False)
        self._launch(
            lambda: compare_scalar_fields(left, right, field).to_dict(),
            self.compare_text,
            lambda value: self.difference_button.setEnabled(value["cells"] <= 250_000),
        )

    def compare_conservative(self):
        left, right, field = (
            self.left_file.text(),
            self.right_file.text(),
            self.compare_field.text(),
        )
        self.difference_button.setEnabled(False)
        self._launch(
            lambda: conservative_box_remap(left, right, field)[1].to_dict(),
            self.compare_text,
        )

    def compare_convex(self):
        left, right, field = (
            self.left_file.text(),
            self.right_file.text(),
            self.compare_field.text(),
        )
        self.difference_button.setEnabled(False)
        self._launch(
            lambda: conservative_convex_remap(left, right, field)[1].to_dict(),
            self.compare_text,
        )

    def request_volume_partition(self):
        path = self.volume_file.text().strip()
        if not path:
            self.volume_note.setText(
                "Choose an HDF5 file containing the native CFDX volume mesh"
            )
            return
        field = self.volume_scalar.text().strip() or None
        self.volume_note.setText(
            "Reading selected volume tile in background; no full mesh allocation"
        )
        self.volume_partition_requested.emit(
            path, self.volume_start.value(), self.volume_count.value(), field
        )

    def compare_telemetry(self):
        left, right = self.reference_run.text(), self.candidate_run.text()
        self._launch(lambda: compare_run_reports(left, right), self.telemetry_text)
