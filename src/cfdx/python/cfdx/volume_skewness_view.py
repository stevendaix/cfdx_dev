"""Interactive cut-plane viewer for volume-cell skewness diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import (
    QComboBox, QDialog, QLabel, QVBoxLayout, QHBoxLayout, QPushButton,
)

from .volume_mesh_quality import VolumeSkewnessReport, analyze_cfdx_volume_skewness


class VolumeSkewnessDialog(QDialog):
    """Show bad cells in red, the cut context transparently, and local neighbours."""

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
            raise RuntimeError("PyVista and PyVistaQt are required for volume inspection") from exc
        self._pv = pv
        self.report: VolumeSkewnessReport = analyze_cfdx_volume_skewness(
            path, expected_signature=expected_signature
        )
        source = Path(path)
        stat = source.stat()
        signature = (int(stat.st_dev), int(stat.st_ino), int(stat.st_size), int(stat.st_mtime_ns))
        if expected_signature is not None and signature != expected_signature:
            raise RuntimeError("mesh file changed before volume rendering")
        with h5py.File(source, "r") as h5:
            points = np.asarray(h5["points"][:], dtype=float)
            face_vertices = np.asarray(h5["face_vertices"][:], dtype=np.int64)
            face_offsets = np.asarray(h5["face_offsets"][:], dtype=np.int64)
            cell_faces = np.asarray(h5["cell_faces"][:], dtype=np.int64)
            cell_offsets = np.asarray(h5["cell_offsets"][:], dtype=np.int64)
        stat_after = source.stat()
        after = (int(stat_after.st_dev), int(stat_after.st_ino), int(stat_after.st_size), int(stat_after.st_mtime_ns))
        if after != signature:
            raise RuntimeError("mesh file changed while volume rendering topology was loaded")
        vtk_cells: list[int] = []
        for cell_id in range(len(cell_offsets) - 1):
            face_ids = cell_faces[cell_offsets[cell_id]:cell_offsets[cell_id + 1]]
            vtk_cells.append(len(face_ids))
            for face_id in face_ids:
                vertices = face_vertices[face_offsets[face_id]:face_offsets[face_id + 1]]
                vtk_cells.append(len(vertices))
                vtk_cells.extend(int(value) for value in vertices)
        grid = pv.UnstructuredGrid(
            np.asarray(vtk_cells, dtype=np.int64),
            np.full(self.report.cell_skewness.size, pv.CellType.POLYHEDRON, dtype=np.uint8),
            points,
        )
        grid.cell_data["skewness"] = self.report.cell_skewness
        grid.cell_data["bad"] = (
            self.report.cell_skewness > self.report.threshold
        ).astype(np.uint8)
        self._grid = grid
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
        self._slice_axis.addItems(["X", "Y", "Z"])
        self._slice_axis.setCurrentText("Z")
        self._slice_button = QPushButton("Update cut", self)
        self._slice_button.clicked.connect(self._render_cut)
        controls = QHBoxLayout()
        controls.addWidget(QLabel("Problem cell:"))
        controls.addWidget(self._cell_selector, 1)
        controls.addWidget(QLabel("Cut normal:"))
        controls.addWidget(self._slice_axis)
        controls.addWidget(self._slice_button)
        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(self._plotter, 1)
        layout.addWidget(self._details)
        self._cell_selector.currentIndexChanged.connect(self._show_selected_cell)
        self._render_cut()
        self._show_selected_cell()

    def _render_cut(self) -> None:
        axis = self._slice_axis.currentText().lower()
        bounds = self._grid.bounds
        index = {"x": 0, "y": 1, "z": 2}[axis]
        origin = list(self._grid.center)
        origin[index] = 0.5 * (bounds[2 * index] + bounds[2 * index + 1])
        normal = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}[axis]
        cut = self._grid.slice(normal=normal, origin=origin)
        self._plotter.clear()
        if cut.n_cells:
            bad_cut = cut.threshold(0.5, scalars="bad")
            good_cut = cut.threshold(0.5, scalars="bad", invert=True)
            if good_cut.n_cells:
                self._plotter.add_mesh(
                    good_cut, color="lightgray", opacity=0.10,
                    show_edges=True, line_width=0.5,
                )
            if bad_cut.n_cells:
                self._plotter.add_mesh(
                    bad_cut, color="red", opacity=1.0, show_edges=True,
                    line_width=1.5,
                )
        self._plotter.add_text(
            f"Cut {axis.upper()} · red: skewness > {self.report.threshold:g} · "
            "transparent: remaining cells",
            font_size=10,
        )
        self._plotter.reset_camera()
        self._plotter.render()
        self._show_selected_cell()

    def _show_selected_cell(self, *_args) -> None:
        cell_id = self._cell_selector.currentData()
        if cell_id is None:
            return
        neighbours = self.report.cell_neighbours[int(cell_id)]
        self._details.setText(
            f"Cell {cell_id} · skewness {self.report.cell_skewness[cell_id]:.6g} · "
            f"threshold {self.report.threshold:g}\n"
            f"Adjacent cells ({len(neighbours)}): "
            + (", ".join(str(value) for value in neighbours) or "none")
        )
        ids = [int(cell_id), *neighbours]
        local = self._grid.extract_cells(ids).extract_surface()
        self._plotter.add_mesh(
            local, color="yellow", opacity=0.95, show_edges=True,
            line_width=2.0, name="selected-cell-neighbours",
        )
        self._plotter.render()

    def closeEvent(self, event) -> None:
        self._plotter.close()
        super().closeEvent(event)
