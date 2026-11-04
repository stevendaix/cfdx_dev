"""Optional PyVistaQt 3D view for the CFDX GUI."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .case_inspector import CaseFileSummary, inspect_cfdx_file
from .mesh_topology import validate_cell_topology, validate_surface_topology
from .selection import parse_selection_id, selection_label
from .visualization_pipeline import VisualizationOperation, VisualizationPipeline

try:
    from PySide6.QtWidgets import QVBoxLayout, QWidget
except ImportError:
    QWidget = None


def _require_previewable_mesh(
    path: Path, summary: CaseFileSummary | None
) -> CaseFileSummary:
    """Fail closed before any eager HDF5 topology arrays are materialized."""
    inspected = summary or inspect_cfdx_file(path)
    if Path(inspected.path).resolve() != path.resolve():
        raise ValueError("mesh preview summary belongs to a different file")
    if inspected.file_signature is None:
        raise ValueError("mesh preview summary has no verified file identity")
    stat = path.stat()
    current_signature = (
        int(stat.st_dev),
        int(stat.st_ino),
        int(stat.st_size),
        int(stat.st_mtime_ns),
    )
    if current_signature != inspected.file_signature:
        raise ValueError(
            "mesh file changed since metadata inspection; inspect it again"
        )
    if not inspected.is_cfdx or not inspected.has_mesh:
        raise ValueError("file does not contain a recognised CFDX mesh")
    if not inspected.preview_supported:
        raise ValueError(
            "full in-memory mesh preview deferred: "
            f"{inspected.preview_block_reason or 'mesh exceeds renderer limits'}"
        )
    return inspected


def _require_previewable_dat_mesh(
    path: Path, summary: CaseFileSummary | None
) -> CaseFileSummary:
    """Preflight cell-to-face arrays before reading a DAT preview topology."""
    inspected = _require_previewable_mesh(path, summary)
    import h5py

    with h5py.File(path, "r") as h5:
        datasets = {}
        for name in ("cell_faces", "cell_offsets"):
            link = h5.get(name, getlink=True)
            dataset = h5.get(name)
            if not isinstance(link, h5py.HardLink) or not isinstance(
                dataset, h5py.Dataset
            ):
                raise TypeError(f"DAT preview requires an internal HDF5 {name} dataset)
            if dataset.is_virtual or dataset.external:
                raise ValueError(f"DAT preview does not support external {name}")
            if dataset.ndim != 1 or dataset.dtype.kind not in "iu":
                raise ValueError(f"DAT preview requires one-dimensional integer {name}")
            datasets[name] = dataset

        if datasets["cell_offsets"].shape != (int(inspected.n_cells) + 1,):
            raise ValueError("DAT preview cell_offsets do not match the cell count")
        connectivity_entries = int(datasets["cell_faces"].shape[0])
        source_bytes = sum(
            int(dataset.size * dataset.dtype.itemsize) for dataset in datasets.values()
        )
        if connectivity_entries > 2_000_000:
            raise ValueError(
                "DAT cell connectivity exceeds the 2M-entry eager-load budget"
            )
        if source_bytes > 32 * 1024 * 1024:
            raise ValueError("DAT cell topology exceeds the 32 MiB eager-load budget")
    return inspected


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
            raise RuntimeError(
                "PyVista and pyvistaqt are required for the 3D GUI"
            ) from exc
        self._pv = pv
        self.plotter = QtInteractor(self)
        self._actors: dict[str, object] = {}
        self._current_result: Path | None = None
        self._selected_field: str | None = None
        self._dat_dataset = None
        self._dat_fields: list[str] = []
        self._visualization_pipeline = VisualizationPipeline()
        self._reset_selection_state()
        layout = QVBoxLayout(self)
        layout.addWidget(self.plotter)

    def _dataset(self):
        if self._dat_dataset is not None:
            return self._dat_dataset
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
        self._visualization_pipeline = VisualizationPipeline()
        self._reset_selection_state()
        self._current_result = path
        self._dat_dataset = None
        self._dat_fields = []
        dataset = self._dataset()
        self.plotter.add_mesh(
            dataset,
            scalars=self._selected_field
            if self._selected_field in self._available_fields(dataset)
            else None,
        )
        self.plotter.reset_camera()
        self.plotter.render()

    @staticmethod
    def patch_actor_id(stable_id: str) -> str:
        """Return the renderer actor ID for a stable application patch ID."""
        if not stable_id.startswith("patch:") or not stable_id.split(":", 1)[1]:
            raise ValueError(f"invalid patch stable ID: {stable_id!r}")
        return stable_id

    def load_cfdx_mesh(
        self, source: str, *, summary: CaseFileSummary | None = None
    ) -> None:
        """Load a guarded CFDX HDF5 mesh as selectable 3D actors."""
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(source)
        inspected = _require_previewable_mesh(path, summary)
        import h5py
        import numpy as np

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
            validate_surface_topology(
                points,
                face_vertices,
                face_offsets,
                patch_face_ids,
                patch_offsets,
                n_faces=int(inspected.n_faces),
                n_patches=len(metadata),
            )
            self.plotter.clear()
            self._actors.clear()
            self._visualization_pipeline = VisualizationPipeline()
            self._reset_selection_state()
            self._mesh_points = points
            for index, item in enumerate(metadata):
                fields = item.split(":")
                if len(fields) < 4:
                    raise ValueError(f"invalid boundary patch metadata: {item!r}")
                name = fields[0]
                start, end = int(patch_offsets[index]), int(patch_offsets[index + 1])
                faces = []
                patch_faces: list[tuple[int, ...]] = []
                for face_id in patch_face_ids[start:end]:
                    fid = int(face_id)
                    begin, stop = int(face_offsets[fid]), int(face_offsets[fid + 1])
                    vertices = tuple(int(v) for v in face_vertices[begin:stop])
                    patch_faces.append(vertices)
                    faces.extend([len(vertices), *vertices])
                actor_id = self.patch_actor_id(f"patch:{name}")
                mesh = self._pv.PolyData(points, np.asarray(faces, dtype=np.int64))
                self._actors[actor_id] = self.plotter.add_mesh(mesh, name=actor_id)
                self._patch_face_vertices[actor_id] = tuple(patch_faces)
                self._actor_regions[actor_id] = "default"
            self.plotter.reset_camera()
            self.plotter.render()

    def load_cfdx_dat(
        self,
        case_path: str,
        dat_path: str,
        *,
        summary: CaseFileSummary | None = None,
    ) -> list[str]:
        """Load a guarded CFDX mesh and DAT checkpoint into one dataset."""
        mesh_path = Path(case_path)
        if not mesh_path.exists():
            raise FileNotFoundError(case_path)
        inspected = _require_previewable_dat_mesh(mesh_path, summary)
        import h5py
        import numpy as np

        from .dat_io import read_dat_restart

        restart = read_dat_restart(dat_path)
        with h5py.File(mesh_path, "r") as h5:
            for name in ("cell_faces", "cell_offsets"):
                link = h5.get(name, getlink=True)
                if not isinstance(link, h5py.HardLink) or not isinstance(
                    h5.get(name), h5py.Dataset
                ):
                    raise TypeError(
                        f"DAT preview requires an internal HDF5 {name} dataset"
                    )
            points = np.asarray(h5["points"][:], dtype=float)
            face_vertices = np.asarray(h5["face_vertices"][:], dtype=np.int64)
            face_offsets = np.asarray(h5["face_offsets"][:], dtype=np.int64)
            cell_faces = np.asarray(h5["cell_faces"][:], dtype=np.int64)
            cell_offsets = np.asarray(h5["cell_offsets"][:], dtype=np.int64)
            patch_face_ids = np.asarray(h5["patch_face_ids"][:], dtype=np.int64)
            patch_offsets = np.asarray(h5["patch_face_offsets"][:], dtype=np.int64)
            raw = h5.attrs.get("boundary_patches", "")
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            metadata = [item for item in str(raw).split(";") if item]

        validate_surface_topology(
            points,
            face_vertices,
            face_offsets,
            patch_face_ids,
            patch_offsets,
            n_faces=int(inspected.n_faces),
            n_patches=len(metadata),
        )
        mesh_cells = int(inspected.n_cells)
        validate_cell_topology(
            cell_faces,
            cell_offsets,
            face_offsets,
            n_cells=mesh_cells,
            n_faces=int(inspected.n_faces),
        )
        if restart.cells != mesh_cells:
            raise ValueError(
                f"DAT cell count {restart.cells} does not match HDF5 mesh cell count {mesh_cells}"
            )

        cells: list[int] = []
        for cell_id in range(mesh_cells):
            first, last = int(cell_offsets[cell_id]), int(cell_offsets[cell_id + 1])
            faces = cell_faces[first:last]
            cells.append(len(faces))
            for face_id in faces:
                fid = int(face_id)
                begin, end = int(face_offsets[fid]), int(face_offsets[fid + 1])
                vertices = face_vertices[begin:end]
                cells.append(len(vertices))
                cells.extend(int(v) for v in vertices)

        celltypes = np.full(mesh_cells, self._pv.CellType.POLYHEDRON, dtype=np.uint8)
        dataset = self._pv.UnstructuredGrid(
            np.asarray(cells, dtype=np.int64), celltypes, points
        )
        for name, field in restart.fields.items():
            values = np.asarray(field.values, dtype=float).reshape(
                mesh_cells, field.dimension
            )
            dataset.cell_data[name] = values[:, 0] if field.dimension == 1 else values

        self.plotter.clear()
        self._actors.clear()
        self._visualization_pipeline = VisualizationPipeline()
        self._reset_selection_state()
        self._current_result = None
        self._selected_field = None
        self._dat_dataset = dataset
        self._dat_fields = list(restart.fields)
        self.plotter.add_mesh(
            dataset, scalars=self._dat_fields[0] if self._dat_fields else None
        )
        self.plotter.reset_camera()
        self.plotter.render()
        return list(self._dat_fields)

    def _render_visualization(self) -> None:
        """Render the replayed pipeline without replacing its source dataset."""
        source = self._dataset()
        derived = self._visualization_pipeline.evaluate(source)
        self.plotter.clear()
        self._selection_highlight = None
        self._selection_id = None
        self.plotter.add_mesh(
            derived,
            scalars=(
                self._selected_field
                if not self._visualization_pipeline.operations
                and self._selected_field in self._available_fields(derived)
                else None
            ),
        )
        self.plotter.render()

    def set_field(self, field: str) -> None:
        """Display a selected result field on the current source dataset."""
        dataset = self._dataset()
        if field not in self._available_fields(dataset):
            raise KeyError(field)
        self._selected_field = field
        self._visualization_pipeline = VisualizationPipeline()
        self._render_visualization()

    def contour(self, scalars: str, isosurfaces: int = 10) -> None:
        """Append a contour operation to the current visualization pipeline."""
        operation = VisualizationOperation.contour(scalars, isosurfaces)
        if scalars not in self._available_fields(self._dataset()):
            raise KeyError(scalars)
        self._visualization_pipeline = self._visualization_pipeline.append(operation)
        self._render_visualization()

    def slice(self, normal: tuple[float, float, float] = (1.0, 0.0, 0.0)) -> None:
        """Append a slice operation to the current visualization pipeline."""
        operation = VisualizationOperation.slice(normal)
        self._visualization_pipeline = self._visualization_pipeline.append(operation)
        self._render_visualization()

    def glyph(self, scalars: str, factor: float = 1.0) -> None:
        """Append a glyph operation to the current visualization pipeline."""
        operation = VisualizationOperation.glyph(scalars, factor)
        if scalars not in self._available_fields(self._dataset()):
            raise KeyError(scalars)
        self._visualization_pipeline = self._visualization_pipeline.append(operation)
        self._render_visualization()

    def remove_visualization_operation(self, index: int) -> None:
        """Remove one pipeline operation and replay the remaining operations."""
        self._visualization_pipeline = self._visualization_pipeline.remove(index)
        self._render_visualization()

    def replace_visualization_operation(
        self, index: int, operation: VisualizationOperation
    ) -> None:
        """Replace one pipeline operation and replay the updated pipeline."""
        if operation.scalars is not None:
            dataset = self._dataset()
            if operation.scalars not in self._available_fields(dataset):
                raise KeyError(operation.scalars)
        self._visualization_pipeline = self._visualization_pipeline.replace(
            index, operation
        )
        self._render_visualization()

    def move_visualization_operation(self, from_index: int, to_index: int) -> None:
        """Move one pipeline operation and replay the reordered pipeline."""
        self._visualization_pipeline = self._visualization_pipeline.move(
            from_index, to_index
        )
        self._render_visualization()

    def reset_visualization(self) -> None:
        """Remove all derived operations and show the original dataset again."""
        self._visualization_pipeline = self._visualization_pipeline.clear()
        self._render_visualization()

    def select(self, stable_id: str) -> None:
        """Highlight one stable selection target (patch, region or face)."""
        target = parse_selection_id(stable_id)
        self._remove_selection_highlight()
        if target.kind == "patch":
            if stable_id not in self._actors:
                raise KeyError(stable_id)
            for actor_id, actor in self._actors.items():
                actor.GetProperty().SetOpacity(1.0 if actor_id == stable_id else 0.2)
        elif target.kind == "region":
            self._select_region(target)
        else:
            self._select_face(target)
        self._selection_id = stable_id
        self.plotter.render()

    def _select_region(self, target: Any) -> None:
        """Highlight every actor of one region and dim the others."""
        region_actors = [
            actor
            for actor_id, actor in self._actors.items()
            if self._actor_regions.get(actor_id) == target.name
        ]
        if not region_actors:
            raise KeyError(target.stable_id)
        for actor in self._actors.values():
            actor.GetProperty().SetOpacity(0.2)
        for actor in region_actors:
            actor.GetProperty().SetOpacity(1.0)

    def _select_face(self, target: Any) -> None:
        """Add a dedicated highlight actor for one face of a loaded patch."""
        import numpy as np

        patch_id = target.stable_id.rsplit("/face:", 1)[0]
        patch_faces = self._patch_face_vertices.get(patch_id)
        if not patch_faces:
            raise KeyError(target.stable_id)
        if target.face_index is None or target.face_index >= len(patch_faces):
            raise KeyError(target.stable_id)
        vertices = patch_faces[target.face_index]
        face_mesh = np.asarray([len(vertices), *vertices], dtype=np.int64)
        mesh = self._pv.PolyData(self._mesh_points, face_mesh)
        for actor in self._actors.values():
            actor.GetProperty().SetOpacity(0.2)
        self._selection_highlight = self.plotter.add_mesh(
            mesh, color="#ff9f43", name="workbench.selection.face"
        )

    def _remove_selection_highlight(self) -> None:
        """Drop the face-highlight actor without touching the mesh actors."""
        if self._selection_highlight is not None:
            self.plotter.remove_actor("workbench.selection.face")
            self._selection_highlight = None

    def _reset_selection_state(self) -> None:
        """Forget every per-dataset selection tracking structure."""
        self._mesh_points = None
        self._patch_face_vertices = {}
        self._actor_regions = {}
        self._selection_id = None
        self._selection_highlight = None

    def selection_id(self) -> str | None:
        """Return the stable id of the currently highlighted target."""
        return self._selection_id

    def selection_summary(self) -> str:
        """Return a human-readable summary for the status bar."""
        if self._selection_id is None:
            return "No selection"
        return selection_label(parse_selection_id(self._selection_id))

    def clear(self) -> None:
        self.plotter.clear()
        self._actors.clear()
        self._current_result = None
        self._selected_field = None
        self._dat_dataset = None
        self._dat_fields = []
        self._visualization_pipeline = VisualizationPipeline()
        self._reset_selection_state()

    def closeEvent(self, event) -> None:
        self.plotter.close()
        super().closeEvent(event)
