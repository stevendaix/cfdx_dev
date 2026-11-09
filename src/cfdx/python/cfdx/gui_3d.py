"""Optional PyVistaQt 3D view for the CFDX GUI."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .case_inspector import CaseFileSummary, inspect_cfdx_file
from .cell_convergence import sampled_cell_scalar_map
from .mesh_topology import validate_cell_topology, validate_surface_topology
from .selection import parse_selection_id, selection_label
from .surface_streaming import read_surface_preview
from .theme import TOKENS
from .visualization_pipeline import VisualizationOperation, VisualizationPipeline
from .volume_partitions import read_volume_partition

try:
    from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal
    from PySide6.QtWidgets import QVBoxLayout, QWidget
except ImportError:
    QWidget = None
    Signal = None
    QObject = None
    QRunnable = None
    QThreadPool = None


class _SurfaceSignals(QObject if QObject is not None else object):
    done = Signal(int, object, object) if Signal is not None else None


class _SurfaceJob(QRunnable if QRunnable is not None else object):
    def __init__(self, revision: int, path: Path, signature: tuple):
        super().__init__()
        self.revision, self.path, self.signature = revision, path, signature
        self.events = _SurfaceSignals()

    def run(self):
        try:
            result = read_surface_preview(self.path)
            stat = self.path.stat()
            current = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)
            if current != self.signature:
                raise RuntimeError(
                    "mesh changed while constructing progressive surface preview"
                )
            self.events.done.emit(self.revision, result, None)
        except Exception as exc:  # noqa: BLE001
            self.events.done.emit(self.revision, None, str(exc))


class _VolumeJob(QRunnable if QRunnable is not None else object):
    def __init__(
        self, revision: int, path: Path, start_cell: int, count: int, field: str | None
    ):
        super().__init__()
        self.revision, self.path, self.start_cell = revision, path, start_cell
        self.count, self.field = count, field
        self.events = _SurfaceSignals()

    def run(self):
        try:
            result = read_volume_partition(
                self.path, self.start_cell, count=self.count, field=self.field
            )
            self.events.done.emit(self.revision, result, None)
        except Exception as exc:  # noqa: BLE001
            self.events.done.emit(self.revision, None, str(exc))


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
                raise TypeError(f"DAT preview requires an internal HDF5 {name} dataset")
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

    selection_changed = Signal(str) if Signal is not None else None
    progressive_preview_ready = Signal(str) if Signal is not None else None
    progressive_preview_failed = Signal(str) if Signal is not None else None
    volume_partition_ready = Signal(str) if Signal is not None else None
    volume_partition_failed = Signal(str) if Signal is not None else None

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
        self.set_workbench_background(TOKENS["surface.canvas"])
        self._actors: dict[str, object] = {}
        self._current_result: Path | None = None
        self._selected_field: str | None = None
        self._dat_dataset = None
        self._dat_fields: list[str] = []
        self._loaded_mesh_for_hotspots: Path | None = None
        self._hotspot_trace = None
        self._hotspot_frame = None
        self._source_dataset = None
        self._rendered_dataset = None
        self._rendered_scalars: str | None = None
        self._visualization_pipeline = VisualizationPipeline()
        self._preview_revision = 0
        self._preview_jobs = {}
        self._preview_pool = QThreadPool.globalInstance()
        self._reset_selection_state()
        layout = QVBoxLayout(self)
        layout.addWidget(self.plotter)

    def set_workbench_background(self, color: str) -> None:
        """Update the VTK scene background; never alter a field colormap."""
        self.plotter.set_background(color)

    def _dataset(self):
        if self._dat_dataset is not None:
            return self._dat_dataset
        if self._source_dataset is not None:
            return self._source_dataset
        if self._current_result is None:
            raise RuntimeError("no result dataset loaded")
        self._source_dataset = self._pv.read(self._current_result)
        return self._source_dataset

    def _available_fields(self, dataset) -> set[str]:
        return set(dataset.point_data.keys()) | set(dataset.cell_data.keys())

    def load(self, source: str) -> None:
        self._preview_revision += 1
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(source)
        # Decode the source before changing the current scene. Reuse it for
        # subsequent filter replays rather than rereading the whole file.
        dataset = self._pv.read(path)
        scalars = (
            self._selected_field
            if self._selected_field in self._available_fields(dataset)
            else None
        )
        self.plotter.disable_picking()
        self.plotter.clear()
        self._actors.clear()
        self._visualization_pipeline = VisualizationPipeline()
        self._reset_selection_state()
        self._current_result = path
        self._loaded_mesh_for_hotspots = None
        self._hotspot_trace = None
        self._dat_dataset = None
        self._source_dataset = dataset
        self._dat_fields = []
        self.plotter.add_mesh(dataset, scalars=scalars)
        self._rendered_dataset = dataset
        self._rendered_scalars = scalars
        self.plotter.reset_camera()
        self.plotter.render()

    def load_progressive_surface(
        self, source: str, *, summary: CaseFileSummary | None = None
    ) -> None:
        """Schedule a safe, decimated surface-only preview off the Qt UI thread.

        This is never a full-volume field load: sample counts/coverage are
        displayed explicitly and the source cell IDs are not manufactured.
        """
        path = Path(source)
        inspected = summary or inspect_cfdx_file(path)
        if not inspected.has_mesh or not inspected.is_cfdx:
            raise ValueError("no recognised CFDX mesh to preview")
        if inspected.file_signature is None:
            raise ValueError("missing mesh identity metadata")
        stat = path.stat()
        current = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)
        if current != inspected.file_signature:
            raise ValueError("mesh changed since metadata inspection")
        self._preview_revision += 1
        revision = self._preview_revision
        job = _SurfaceJob(revision, path, current)
        self._preview_jobs[revision] = job
        job.events.done.connect(self._on_surface_preview)
        self._preview_pool.start(job)

    def load_progressive_volume(
        self,
        source: str,
        *,
        start_cell: int = 0,
        count: int = 128,
        field: str | None = None,
    ) -> None:
        """Schedule one true volume tile, never materialize the entire mesh."""
        path = Path(source)
        if not path.is_file() or path.is_symlink():
            raise ValueError("volume tile source must be a regular local HDF5 file")
        self._preview_revision += 1
        revision = self._preview_revision
        job = _VolumeJob(revision, path, start_cell, count, field)
        self._preview_jobs[revision] = job
        job.events.done.connect(self._on_volume_partition)
        self._preview_pool.start(job)

    def _on_volume_partition(self, revision: int, piece, error) -> None:
        self._preview_jobs.pop(revision, None)
        if revision != self._preview_revision:
            return  # Newer file/project/partition invalidates the old tile.
        if error:
            self.volume_partition_failed.emit(str(error))
            return
        try:
            from pathlib import Path

            stat = Path(piece.source).stat()
            if (
                stat.st_dev,
                stat.st_ino,
                stat.st_size,
                stat.st_mtime_ns,
            ) != piece.source_identity:
                raise RuntimeError("volume tile source changed before rendering")
            tile = piece.as_pyvista(self._pv)
            self.plotter.disable_picking()
            self.plotter.clear()
            self._actors.clear()
            self._reset_selection_state()
            self._hotspot_trace = None
            self._loaded_mesh_for_hotspots = None  # Partial cells != all cell IDs.
            self._current_result = None
            self._dat_dataset = None
            self._source_dataset = tile
            self._rendered_dataset = tile
            self._dat_fields = [piece.field_name] if piece.field_name else []
            self._selected_field = piece.field_name
            self._rendered_scalars = piece.field_name
            self._visualization_pipeline = VisualizationPipeline()
            self.plotter.add_mesh(
                tile,
                scalars=piece.field_name,
                name="workbench.volume_partition",
                show_edges=True,
            )
            self.plotter.reset_camera()
            self.plotter.render()
        except Exception as exc:  # noqa: BLE001
            self.volume_partition_failed.emit(str(exc))
            return
        first, last = piece.original_cell_ids[[0, -1]]
        self.volume_partition_ready.emit(
            f"Real volume tile {first}-{last} of {piece.total_cells:,} CFDX cells "
            f"({piece.fraction:.2%} of source, not full-volume or GPU qualification)"
        )

    def _on_surface_preview(self, revision: int, preview, error) -> None:
        self._preview_jobs.pop(revision, None)
        if revision != self._preview_revision:
            return  # A newer project/file was opened; never replace its scene.
        if error:
            self.progressive_preview_failed.emit(str(error))
            return
        try:
            # Prepare VTK objects before clearing the previously valid scene.
            mesh = self._pv.PolyData(preview.points, preview.faces)
            mesh.cell_data["CFDX_sampled_patch_index"] = preview.patch_indices
            self.plotter.disable_picking()
            self.plotter.clear()
            self._actors.clear()
            self._reset_selection_state()
            self._loaded_mesh_for_hotspots = None
            self._hotspot_trace = None
            self._current_result = None
            self._selected_field = None
            self._dat_dataset = None
            self._source_dataset = None  # Sampling is not a full analysis field.
            self._dat_fields = []
            self._visualization_pipeline = VisualizationPipeline()
            self._rendered_dataset = mesh
            self._rendered_scalars = "CFDX_sampled_patch_index"
            self.plotter.add_mesh(
                mesh,
                scalars="CFDX_sampled_patch_index",
                cmap="tab20",
                name="workbench.progressive_surface_sample",
                show_scalar_bar=True,
                scalar_bar_args={"title": "Sampled patch ID"},
            )
            self.plotter.reset_camera()
            self.plotter.render()
            self.progressive_preview_ready.emit(
                f"Surface LOD: {len(preview.original_face_ids):,} / "
                f"{preview.n_boundary_faces:,} boundary faces "
                f"({preview.coverage:.2%}); NO full-volume fields or cell map"
            )
        except Exception as exc:  # noqa: BLE001
            self.progressive_preview_failed.emit(
                f"surface preview rendering failed: {exc}"
            )

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
        if self is not None:
            self._preview_revision += 1
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
            self.plotter.disable_picking()
            self.plotter.clear()
            self._actors.clear()
            self._visualization_pipeline = VisualizationPipeline()
            self._source_dataset = None
            self._rendered_dataset = None
            self._rendered_scalars = None
            self._dat_dataset = None
            self._current_result = None
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
            self._loaded_mesh_for_hotspots = path.resolve()
            self._hotspot_trace = None
            self.plotter.reset_camera()
            self.plotter.enable_mesh_picking(
                callback=self._on_mesh_picked,
                use_actor=True,
                show=False,
            )
            self.plotter.render()

    def _stable_id_for_actor(self, picked_actor: object) -> str | None:
        """Resolve a picked VTK/PyVista actor to its registered mesh patch ID."""
        name = getattr(picked_actor, "name", None)
        if isinstance(name, str) and name in self._actors:
            return name
        for stable_id, actor in self._actors.items():
            if picked_actor is actor:
                return stable_id
            try:
                if picked_actor == actor:
                    return stable_id
            except (TypeError, ValueError):
                continue
        return None

    def _on_mesh_picked(self, picked_actor: object) -> None:
        """Publish patch picks; never guess a face/cell ID from VTK order."""
        stable_id = self._stable_id_for_actor(picked_actor)
        if stable_id is not None and self.selection_changed is not None:
            self.selection_changed.emit(stable_id)

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
        if self is not None:
            self._preview_revision += 1
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

        self.plotter.disable_picking()
        self.plotter.clear()
        self._actors.clear()
        self._visualization_pipeline = VisualizationPipeline()
        self._reset_selection_state()
        self._current_result = None
        self._selected_field = None
        self._dat_dataset = dataset
        self._loaded_mesh_for_hotspots = mesh_path.resolve()
        self._hotspot_trace = None
        self._source_dataset = None
        self._dat_fields = list(restart.fields)
        scalars = self._dat_fields[0] if self._dat_fields else None
        self.plotter.add_mesh(dataset, scalars=scalars)
        self._rendered_dataset = dataset
        self._rendered_scalars = scalars
        self.plotter.reset_camera()
        self.plotter.render()
        return list(self._dat_fields)

    def _render_visualization(
        self, candidate: VisualizationPipeline | None = None
    ) -> None:
        """Stage derived data first; commit filter state only after rendering.

        If a filter rejects its input, the existing VTK scene is untouched.
        If VTK rendering fails, try restoring the last successfully displayed
        dataset, without ever publishing the failed operation to the pipeline.
        """
        pipeline = candidate if candidate is not None else self._visualization_pipeline
        source = self._dataset()
        derived = pipeline.evaluate(source)
        scalars = (
            self._selected_field
            if not pipeline.operations
            and self._selected_field in self._available_fields(derived)
            else None
        )
        old_display = self._rendered_dataset
        old_scalars = self._rendered_scalars
        self.plotter.disable_picking()
        try:
            self.plotter.clear()
            self.plotter.add_mesh(derived, scalars=scalars)
            self.plotter.render()
        except Exception:  # noqa: BLE001, S110
            if old_display is not None:
                try:
                    self.plotter.clear()
                    self.plotter.add_mesh(old_display, scalars=old_scalars)
                    self.plotter.render()
                except Exception:
                    pass  # Renderer failure still propagates; never claim success.
            raise
        self._selection_highlight = None
        self._selection_id = None
        self._visualization_pipeline = pipeline
        self._rendered_dataset = derived
        self._rendered_scalars = scalars

    def _hotspot_overlay_dataset(self, trace):
        """Return only the unfiltered DAT volume retaining original cell IDs.

        Matching ``n_cells`` is not sufficient for a filtered or re-ordered VTK
        dataset: the scalar array would otherwise color the wrong CFD cells.
        """
        source = self._dat_dataset
        if source is not None and source.n_cells == trace.n_cells:
            return source
        return None

    def show_cell_hotspots(self, trace, frame, metric: str = "momentum") -> None:
        """Overlay a blue→red sampled convergence map, or fall back to centers.

        The UI never invents values for unsampled cells: any missing cell in the
        bounded top-K snapshot remains neutral grey on a cell map, or absent from
        the fallback point view.
        """
        import numpy as np

        mesh_path = self._loaded_mesh_for_hotspots
        if mesh_path is None:
            raise RuntimeError(
                "Load the source CFDX mesh or matching DAT checkpoint before locating hotspots"
            )
        if mesh_path != Path(trace.mesh_source).resolve():
            raise ValueError(
                "Hotspot source mesh differs from the currently displayed mesh"
            )
        from .cell_convergence import cell_trace_metric_name

        title = cell_trace_metric_name(metric)
        samples = frame.ranked(metric, 72)
        self.plotter.remove_actor("workbench.cell_hotspots", render=False)
        self.plotter.remove_actor("workbench.cell_hotspots.map", render=False)
        if not samples:
            self.plotter.render()
            return
        overlay = self._hotspot_overlay_dataset(trace)
        # A large unstructured grid can be many gigabytes. Never make a deep
        # clone merely to add a visualization array: the *read-only* geometry
        # is shared by a shallow copy, and a hard budget protects cell data.
        if overlay is not None and trace.n_cells <= 1_000_000:
            _, sparse_map = sampled_cell_scalar_map(trace, frame, metric)
            scalars = np.asarray(sparse_map, dtype=float)
            # Keep a diagnostic array distinct from the native CFD fields.
            overlay_copy = overlay.copy(deep=False)
            field_name = f"CFDX hotspot {metric}"
            overlay_copy.cell_data[field_name] = scalars
            finite = scalars[np.isfinite(scalars)]
            clim = None
            if finite.size:
                lower, upper = float(np.min(finite)), float(np.max(finite))
                if upper <= lower:
                    # At 1e-10, expanding to 1.0 would hide every variation.
                    # Widen at the observed numerical scale instead.
                    margin = max(abs(lower) * 0.05, 1e-15)
                    lower = max(0.0, lower - margin)
                    upper += margin
                clim = (lower, upper)
            self.plotter.add_mesh(
                overlay_copy,
                scalars=field_name,
                name="workbench.cell_hotspots.map",
                cmap="coolwarm",
                show_scalar_bar=True,
                reset_camera=False,
                scalar_bar_args={"title": title + " (sampled only)"},
                nan_color="#5b6472",
                opacity=0.92,
                clim=clim,
            )
            centers = np.asarray([cell.center for cell in samples], dtype=float)
            points = self._pv.PolyData(centers)
            points.point_data["CFDX_cell_id"] = np.asarray(
                [cell.cell_id for cell in samples], dtype=np.int64
            )
            self.plotter.add_mesh(
                points,
                name="workbench.cell_hotspots",
                color="#ffffff",
                render_points_as_spheres=True,
                point_size=6,
                reset_camera=False,
            )
        else:
            centers = np.asarray([cell.center for cell in samples], dtype=float)
            points = self._pv.PolyData(centers)
            field = {
                "momentum": "momentum_scaled",
                "mass": "mass_imbalance_abs",
                "persistence": "nondecrease_streak",
            }[metric]
            points.point_data[field] = np.asarray(
                [getattr(cell, field) for cell in samples], dtype=float
            )
            points.point_data["CFDX_cell_id"] = np.asarray(
                [cell.cell_id for cell in samples], dtype=np.int64
            )
            self.plotter.add_mesh(
                points,
                scalars=field,
                name="workbench.cell_hotspots",
                cmap="coolwarm",
                render_points_as_spheres=True,
                point_size=14,
                show_scalar_bar=True,
                scalar_bar_args={"title": title + " (sampled only)"},
                reset_camera=False,
            )
        self._hotspot_trace = trace
        self._hotspot_frame = frame
        self.plotter.render()

    def highlight_cell_hotspot(self, trace, cell) -> None:
        """Highlight an original solver cell ID without guessing a VTK index."""
        import numpy as np

        mesh_path = self._loaded_mesh_for_hotspots
        if mesh_path is None or mesh_path != Path(trace.mesh_source).resolve():
            raise ValueError(
                "Open the matching original CFDX mesh before selecting the cell"
            )
        center = np.asarray([cell.center], dtype=float)
        point = self._pv.PolyData(center)
        self.plotter.add_mesh(
            point,
            name="workbench.cell_hotspot.selected",
            color="#ffffff",
            render_points_as_spheres=True,
            point_size=24,
            reset_camera=False,
        )
        # Move the camera toward the actual cell position while preserving its
        # viewing direction. This is not a guessed VTK cell index: the overlay
        # uses the source-mesh coordinate exported by the CFDX solver.
        camera = self.plotter.camera
        old_position = np.asarray(camera.position, dtype=float)
        old_focal = np.asarray(camera.focal_point, dtype=float)
        direction = old_position - old_focal
        length = float(np.linalg.norm(direction))
        if length <= 1e-12:
            direction = np.asarray([1.0, 1.0, 1.0])
            length = float(np.linalg.norm(direction))
        bounds = self.plotter.bounds
        span = np.asarray(
            (bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]),
            dtype=float,
        )
        diagonal = float(np.linalg.norm(span))
        distance = max(diagonal * 0.08, 1e-6)
        camera.position = tuple(center[0] + (direction / length) * distance)
        camera.focal_point = tuple(center[0])
        self.plotter.reset_camera_clipping_range()
        self.plotter.render()

    def show_full_cell_map(self, payload) -> None:
        """Display stable, quantitative blue/red scale on original cell IDs only."""
        import numpy as np

        mesh = self._loaded_mesh_for_hotspots
        if mesh is None or mesh != Path(payload["mesh_source"]).resolve():
            raise RuntimeError(
                "load the exact originating CFDX mesh before showing full-cell map"
            )
        native = self._hotspot_overlay_dataset(
            type("_Cells", (), {"n_cells": payload["n_cells"]})()
        )
        if native is None:
            raise RuntimeError(
                "full field colouring requires a compatible volumetric DAT dataset"
            )
        values = np.asarray(payload["values"], dtype=np.float64)
        if values.shape != (native.n_cells,):
            raise ValueError("full-cell map field does not match the source volume")
        overlay = native.copy(deep=False)
        overlay.cell_data["CFDX_full_cell_convergence"] = values
        self.plotter.add_mesh(
            overlay,
            scalars="CFDX_full_cell_convergence",
            name="workbench.cell_hotspots.full_map",
            cmap="coolwarm",
            nan_color="#59616e",
            show_scalar_bar=True,
            scalar_bar_args={"title": payload["title"]},
            clim=payload["clim"],
            reset_camera=False,
        )
        self.plotter.render()

    def show_simulation_difference(self, left: str, right: str, field: str) -> None:
        """Only a small identical-grid field may be rendered as signed delta."""
        import h5py
        import numpy as np

        from .simulation_compare import compare_scalar_fields

        report = compare_scalar_fields(left, right, field)
        if report.cells > 250_000:
            raise RuntimeError(
                "volume difference display limited to 250,000 cells; use streaming statistics"
            )
        if self._loaded_mesh_for_hotspots is None or self._dat_dataset is None:
            raise RuntimeError("load the original matching DAT volumetric mesh first")
        if self._loaded_mesh_for_hotspots.resolve() != Path(left).resolve():
            raise RuntimeError("displayed geometry is not the reference mesh")
        if self._dat_dataset.n_cells != report.cells:
            raise RuntimeError("VTK cell count does not match verified reference field")
        with h5py.File(left, "r") as a, h5py.File(right, "r") as b:
            lhs = np.asarray(a[field][:], dtype=np.float64)
            rhs = np.asarray(b[field][:], dtype=np.float64)
        delta = rhs - lhs
        mesh = self._dat_dataset.copy(deep=False)
        mesh.cell_data["CFDX_signed_difference"] = delta
        self.plotter.add_mesh(
            mesh,
            scalars="CFDX_signed_difference",
            cmap="coolwarm",
            name="workbench.compare.signed_delta",
            scalar_bar_args={"title": f"Δ {report.field} [{report.units}]"},
            reset_camera=False,
        )
        self.plotter.render()

    def set_field(self, field: str) -> None:
        """Display a selected result field on the current source dataset."""
        dataset = self._dataset()
        if field not in self._available_fields(dataset):
            raise KeyError(field)
        previous = self._selected_field
        self._selected_field = field
        try:
            self._render_visualization(VisualizationPipeline())
        except Exception:
            self._selected_field = previous
            raise

    def contour(self, scalars: str, isosurfaces: int = 10) -> None:
        """Append a contour operation to the current visualization pipeline."""
        self._render_visualization(
            self._visualization_pipeline.append(
                VisualizationOperation.contour(scalars, isosurfaces)
            )
        )

    def slice(self, normal: tuple[float, float, float] = (1.0, 0.0, 0.0)) -> None:
        """Append a slice operation to the current visualization pipeline."""
        self._render_visualization(
            self._visualization_pipeline.append(VisualizationOperation.slice(normal))
        )

    def glyph(self, scalars: str, factor: float = 1.0) -> None:
        """Append a glyph operation to the current visualization pipeline."""
        self._render_visualization(
            self._visualization_pipeline.append(
                VisualizationOperation.glyph(scalars, factor)
            )
        )

    def clip(
        self,
        normal: tuple[float, float, float] = (1.0, 0.0, 0.0),
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
        *,
        invert: bool = True,
    ) -> None:
        self._render_visualization(
            self._visualization_pipeline.append(
                VisualizationOperation.clip(normal, origin, invert=invert)
            )
        )

    def threshold(
        self, scalars: str, minimum: float, maximum: float, *, invert: bool = False
    ) -> None:
        self._render_visualization(
            self._visualization_pipeline.append(
                VisualizationOperation.threshold(
                    scalars, minimum, maximum, invert=invert
                )
            )
        )

    def remove_visualization_operation(self, index: int) -> None:
        """Remove one operation only after the replacement scene is prepared."""
        self._render_visualization(self._visualization_pipeline.remove(index))

    def replace_visualization_operation(
        self, index: int, operation: VisualizationOperation
    ) -> None:
        self._render_visualization(
            self._visualization_pipeline.replace(index, operation)
        )

    def move_visualization_operation(self, from_index: int, to_index: int) -> None:
        self._render_visualization(
            self._visualization_pipeline.move(from_index, to_index)
        )

    def reset_visualization(self) -> None:
        """Remove derived operations; keep the original source dataset."""
        self._render_visualization(self._visualization_pipeline.clear())

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
        self._preview_revision += 1  # Reject in-flight worker results.
        self._loaded_mesh_for_hotspots = None
        self._hotspot_trace = None
        self._hotspot_frame = None
        self.plotter.disable_picking()
        self.plotter.clear()
        self._actors.clear()
        self._current_result = None
        self._selected_field = None
        self._dat_dataset = None
        self._source_dataset = None
        self._rendered_dataset = None
        self._rendered_scalars = None
        self._dat_fields = []
        self._visualization_pipeline = VisualizationPipeline()
        self._reset_selection_state()

    def closeEvent(self, event) -> None:
        self._preview_revision += 1
        self.plotter.close()
        super().closeEvent(event)
