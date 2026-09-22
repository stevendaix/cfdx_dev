"""Optional PyVistaQt 3D view for the CFDX GUI."""
from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    from PySide6.QtWidgets import QVBoxLayout, QWidget
except ImportError:
    QWidget = None


class PyVistaQtView(QWidget if QWidget is not None else object):
    """Qt container around pyvistaqt for CFDX result visualization."""

    def __init__(self, parent: Any = None) -> None:
        if QWidget is None:
            raise RuntimeError("PySide6 is required for the 3D GUI")
        super().__init__(parent)
        try:
            import pyvista as pv
            from pyvistaqt import QtInteractor
        except ImportError as exc:
            raise RuntimeError("PyVista and pyvistaqt are required for the 3D GUI") from exc
        self._pv = pv
        self.plotter = QtInteractor(self)
        self._actors: dict[str, object] = {}
        self._current_result: Path | None = None
        self._selected_field: str | None = None
        layout = QVBoxLayout(self)
        layout.addWidget(self.plotter)

    def _dataset(self):
        if self._current_result is None:
            raise RuntimeError("no result dataset loaded")
        return self._pv.read(self._current_result)

    def _available_fields(self, dataset) -> set[str]:
        return set(dataset.point_data.keys()) | set(dataset.cell_data.keys())

    def load(self, source: str) -> None:
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(source)
        self.plotter.clear()
        self._actors.clear()
        self._current_result = path
        dataset = self._dataset()
        self.plotter.add_mesh(dataset, scalars=self._selected_field if self._selected_field in self._available_fields(dataset) else None)
        self.plotter.reset_camera()
        self.plotter.render()

    def load_cfdx_mesh(self, source: str) -> None:
        """Load real CFDX HDF5 mesh patches as selectable 3D actors."""
        import h5py
        import numpy as np

        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(source)
        with h5py.File(path, "r") as h5:
            points = np.asarray(h5["points"][:], dtype=float)
            face_vertices = np.asarray(h5["face_vertices"][:], dtype=np.int64)
            face_offsets = np.asarray(h5["face_offsets"][:], dtype=np.int64)
            patch_face_ids = np.asarray(h5["patch_face_ids"][:], dtype=np.int64)
            patch_offsets = np.asarray(h5["patch_face_offsets"][:], dtype=np.int64)
            raw = h5.attrs.get("boundary_patches", "")
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            metadata = [item for item in str(raw).split(";") if item]
            if len(patch_offsets) != len(metadata) + 1:
                raise ValueError("boundary patch metadata/offset count mismatch")
            self.plotter.clear()
            self._actors.clear()
            for index, item in enumerate(metadata):
                fields = item.split(":")
                if len(fields) < 4:
                    raise ValueError(f"invalid boundary patch metadata: {item!r}")
                name = fields[0]
                start, end = int(patch_offsets[index]), int(patch_offsets[index + 1])
                faces = []
                for face_id in patch_face_ids[start:end]:
                    fid = int(face_id)
                    begin, stop = int(face_offsets[fid]), int(face_offsets[fid + 1])
                    vertices = face_vertices[begin:stop]
                    faces.extend([len(vertices), *map(int, vertices)])
                actor_id = f"patch:{index}"
                mesh = self._pv.PolyData(points, np.asarray(faces, dtype=np.int64))
                self._actors[actor_id] = self.plotter.add_mesh(mesh, name=actor_id)
            self.plotter.reset_camera()
            self.plotter.render()

    def set_field(self, field: str) -> None:
        """Display a selected result field on the current dataset."""
        dataset = self._dataset()
        if field not in self._available_fields(dataset):
            raise KeyError(field)
        self._selected_field = field
        self.plotter.clear()
        self.plotter.add_mesh(dataset, scalars=field)
        self.plotter.render()

    def contour(self, scalars: str, isosurfaces: int = 10) -> None:
        dataset = self._dataset()
        if scalars not in self._available_fields(dataset):
            raise KeyError(scalars)
        self.plotter.clear()
        self.plotter.add_mesh(dataset.contour(isosurfaces=isosurfaces, scalars=scalars))
        self.plotter.render()

    def slice(self, normal: tuple[float, float, float] = (1.0, 0.0, 0.0)) -> None:
        dataset = self._dataset()
        self.plotter.clear()
        self.plotter.add_mesh(dataset.slice(normal=normal))
        self.plotter.render()

    def glyph(self, scalars: str, factor: float = 1.0) -> None:
        dataset = self._dataset()
        if scalars not in self._available_fields(dataset):
            raise KeyError(scalars)
        array = dataset.point_data.get(scalars)
        if array is None:
            array = dataset.cell_data.get(scalars)
        if array is None:
            raise KeyError(scalars)
        self.plotter.clear()
        kwargs = {"scale": scalars, "factor": factor}
        if getattr(array, "ndim", 1) == 2 and array.shape[1] == 3:
            kwargs["orient"] = scalars
        self.plotter.add_mesh(dataset.glyph(**kwargs))
        self.plotter.render()

    def select(self, object_id: str) -> None:
        """Highlight a stable patch actor by application selection ID."""
        if object_id not in self._actors:
            raise KeyError(object_id)
        for actor_id, actor in self._actors.items():
            actor.GetProperty().SetOpacity(1.0 if actor_id == object_id else 0.2)
        self.plotter.render()

    def clear(self) -> None:
        self.plotter.clear()
        self._actors.clear()
        self._current_result = None
        self._selected_field = None

    def closeEvent(self, event) -> None:
        self.plotter.close()
        super().closeEvent(event)
