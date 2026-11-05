"""Interactive cut-plane viewer for volume-cell skewness diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QVBoxLayout,
)

from .volume_mesh_quality import VolumeSkewnessReport, analyze_cfdx_volume_skewness


class VolumeSkewnessDialog(QDialog):
    """Inspect a mesh cut and highlight the selected cell and its neighbours."""

    def __init__(
        self, path: str | Path, *, expected_signature=None, parent=None
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Volumetric mesh skewness")
        self.resize(1100, 760)
        try:
            import h5py
            import pyvista as pv
            from pyvistaqt import QtInteractor
        except ImportError as exc:
            raise RuntimeError(
                "PyVista and PyVistaQt are required for volume inspection"
            ) from exc
        self._pv = pv
        self.report: VolumeSkewnessReport = analyze_cfdx_volume_skewness(
            path, expected_signature=expected_signature
        )
        source = Path(path)
        stat = source.stat()
        signature = (
            int(stat.st_dev),
            int(stat.st_ino),
            int(stat.st_size),
            int(stat.st_mtime_ns),
        )
        if expected_signature is not None and signature != expected_signature:
            raise RuntimeError("mesh file changed before volume rendering")
        with h5py.File(source, "r") as h5:
            points = np.asarray(h5["points"][:], dtype=float)
            face_vertices = np.asarray(h5["face_vertices"][:], dtype=np.int64)
            face_offsets = np.asarray(h5["face_offsets"][:], dtype=np.int64)
            cell_faces = np.asarray(h5["cell_faces"][:], dtype=np.int64)
            cell_offsets = np.asarray(h5["cell_offsets"][:], dtype=np.int64)
        stat_after = source.stat()
        after = (
            int(stat_after.st_dev),
            int(stat_after.st_ino),
            int(stat_after.st_size),
            int(stat_after.st_mtime_ns),
        )
        if after != signature:
            raise RuntimeError(
                "mesh file changed while volume rendering topology was loaded"
            )
        vtk_cells: list[int] = []
        for cell_id in range(len(cell_offsets) - 1):
            face_ids = cell_faces[cell_offsets[cell_id] : cell_offsets[cell_id + 1]]
            vtk_cells.append(len(face_ids))
            for face_id in face_ids:
                vertices = face_vertices[
                    face_offsets[face_id] : face_offsets[face_id + 1]
                ]
                vtk_cells.append(len(vertices))
                vtk_cells.extend(int(value) for value in vertices)
        grid = pv.UnstructuredGrid(
            np.asarray(vtk_cells, dtype=np.int64),
            np.full(
                self.report.cell_skewness.size, pv.CellType.POLYHEDRON, dtype=np.uint8
            ),
            points,
        )
        grid.cell_data["skewness"] = self.report.cell_skewness
        grid.cell_data["bad"] = (
            self.report.cell_skewness > self.report.threshold
        ).astype(np.uint8)
        self._grid = grid
        self._threshold = float(self.report.threshold)
        self._cut_axis = "z"
        self._cut_normal = (0.0, 0.0, 1.0)
        self._cut_origin = np.asarray(grid.center, dtype=float)
        self._camera_initialized = False
        self._plotter = QtInteractor(self)
        self._cell_selector = QComboBox(self)
        self._cell_selector.setObjectName("volume_skewness.cell_selector")
        for cell_id in self.report.worst_cell_ids:
            self._cell_selector.addItem(
                f"Cell {cell_id} · skewness {self.report.cell_skewness[cell_id]:.4g}",
                cell_id,
            )
        self._details = QLabel(self)
        self._details.setWordWrap(True)
        self._slice_axis = QComboBox(self)
        self._slice_axis.setObjectName("volume_skewness.slice_axis")
        self._slice_axis.addItems(["X", "Y", "Z"])
        self._slice_axis.setCurrentText("Z")
        self._slice_position = QSlider(Qt.Orientation.Horizontal, self)
        self._slice_position.setObjectName("volume_skewness.slice_position")
        self._slice_position.setRange(0, 1000)
        self._slice_position.setValue(500)
        self._slice_position.setToolTip("Move the cut plane through the mesh")
        self._slice_position_label = QLabel(self)
        self._threshold_spin = QDoubleSpinBox(self)
        self._threshold_spin.setObjectName("volume_skewness.threshold")
        self._threshold_spin.setDecimals(3)
        self._threshold_spin.setRange(0.001, 10.0)
        self._threshold_spin.setSingleStep(0.05)
        self._threshold_spin.setValue(self._threshold)
        self._threshold_spin.setToolTip("Cells above this skewness are marked red")
        self._slice_button = QPushButton("Reset cut", self)
        self._slice_button.clicked.connect(self._reset_cut)
        controls = QHBoxLayout()
        controls.addWidget(QLabel("Worst cell:"))
        controls.addWidget(self._cell_selector, 1)
        controls.addWidget(QLabel("Cut normal:"))
        controls.addWidget(self._slice_axis)
        controls.addWidget(QLabel("Position:"))
        controls.addWidget(self._slice_position, 2)
        controls.addWidget(self._slice_position_label)
        controls.addWidget(QLabel("Skewness limit:"))
        controls.addWidget(self._threshold_spin)
        controls.addWidget(self._slice_button)
        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(self._plotter, 1)
        layout.addWidget(self._details)
        self._cell_selector.currentIndexChanged.connect(self._show_selected_cell)
        self._slice_axis.currentTextChanged.connect(self._render_cut)
        self._slice_position.valueChanged.connect(self._render_cut)
        self._threshold_spin.valueChanged.connect(self._update_threshold)
        self._render_cut()

    def _reset_cut(self) -> None:
        self._slice_position.setValue(500)

    def _update_threshold(self, value: float) -> None:
        self._threshold = float(value)
        self._grid.cell_data["bad"] = (
            self.report.cell_skewness > self._threshold
        ).astype(np.uint8)
        self._render_cut()

    def _render_cut(self, *_args) -> None:
        axis = self._slice_axis.currentText().lower()
        bounds = self._grid.bounds
        index = {"x": 0, "y": 1, "z": 2}[axis]
        minimum, maximum = bounds[2 * index], bounds[2 * index + 1]
        fraction = self._slice_position.value() / 1000.0
        position = minimum + fraction * (maximum - minimum)
        origin = np.asarray(self._grid.center, dtype=float)
        origin[index] = position
        normal = {"x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0), "z": (0.0, 0.0, 1.0)}[
            axis
        ]
        self._cut_axis = axis
        self._cut_normal = normal
        self._cut_origin = origin
        self._slice_position_label.setText(f"{position:.5g}")
        cut = self._grid.slice(normal=normal, origin=origin)
        self._plotter.clear()
        if cut.n_cells:
            bad_cut = cut.threshold(0.5, scalars="bad")
            good_cut = cut.threshold(0.5, scalars="bad", invert=True)
            if good_cut.n_cells:
                self._plotter.add_mesh(
                    good_cut,
                    color="lightgray",
                    opacity=0.10,
                    show_edges=True,
                    line_width=0.5,
                    name="volume-cut-good",
                )
            if bad_cut.n_cells:
                self._plotter.add_mesh(
                    bad_cut,
                    color="red",
                    opacity=1.0,
                    show_edges=True,
                    line_width=1.5,
                    name="volume-cut-bad",
                )
        self._plotter.add_text(
            f"Cut {axis.upper()} · red: skewness > {self._threshold:g} · "
            "transparent: remaining cells",
            font_size=10,
        )
        self._show_selected_cell()
        if not self._camera_initialized:
            self._plotter.reset_camera()
            self._camera_initialized = True
        self._plotter.render()

    def _show_selected_cell(self, *_args) -> None:
        cell_id = self._cell_selector.currentData()
        if cell_id is None:
            self._details.setText("No internal-face cell is available for inspection.")
            return
        cell_id = int(cell_id)
        neighbours = self.report.cell_neighbours[cell_id]
        selected_score = self.report.cell_skewness[cell_id]
        selected_non_orth = self.report.cell_non_orthogonality_deg[cell_id]
        self._details.setText(
            f"Cell {cell_id} · skewness {selected_score:.6g} · "
            f"non-orthogonality {selected_non_orth:.3f} deg · "
            f"threshold {self.report.threshold:g}\n"
            f"Adjacent cells ({len(neighbours)}): "
            + (", ".join(str(value) for value in neighbours) or "none")
        )
        selected_cut = self._grid.extract_cells([cell_id]).slice(
            normal=self._cut_normal, origin=self._cut_origin
        )
        if selected_cut.n_cells:
            self._plotter.add_mesh(
                selected_cut,
                color="orange",
                opacity=1.0,
                show_edges=True,
                line_width=3.0,
                name="selected-cell-cut",
            )
        if neighbours:
            neighbours_cut = self._grid.extract_cells(list(neighbours)).slice(
                normal=self._cut_normal, origin=self._cut_origin
            )
            if neighbours_cut.n_cells:
                self._plotter.add_mesh(
                    neighbours_cut,
                    color="yellow",
                    opacity=1.0,
                    show_edges=True,
                    line_width=2.0,
                    name="selected-cell-neighbours-cut",
                )
        if not selected_cut.n_cells:
            self._details.setText(
                self._details.text()
                + "\nSelected cell does not intersect the current cut plane."
            )

    def closeEvent(self, event) -> None:
        self._plotter.close()
        super().closeEvent(event)
