"""Bounded volumetric mesh skewness diagnostics for CFDX HDF5 meshes.

The metric is a face-based internal-face skewness: distance between the face
centroid and the owner/neighbour-centroid line, normalized by the face's
equivalent-area length. Cell scores are the maximum incident internal-face
score. Boundary faces are excluded because they have no neighbour centroid.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import prod
from pathlib import Path

import numpy as np

MAX_VOLUME_CELLS = 100_000
MAX_VOLUME_FACES = 500_000
MAX_VOLUME_CONNECTIVITY = 2_000_000
MAX_VOLUME_SOURCE_BYTES = 64 * 1024 * 1024


@dataclass(frozen=True)
class VolumeSkewnessReport:
    """Per-cell skewness, internal-face scores and cell adjacency."""

    cell_skewness: np.ndarray
    face_skewness: np.ndarray
    cell_neighbours: tuple[tuple[int, ...], ...]
    worst_cell_ids: tuple[int, ...]
    bad_cell_ids: tuple[int, ...]
    threshold: float
    internal_face_count: int
    max_skewness: float
    p95_skewness: float


def analyze_volume_skewness(
    points: np.ndarray,
    face_vertices: np.ndarray,
    face_offsets: np.ndarray,
    owner: np.ndarray,
    neighbour: np.ndarray,
    cell_faces: np.ndarray,
    cell_offsets: np.ndarray,
    *,
    threshold: float = 0.5,
    max_cells: int = MAX_VOLUME_CELLS,
    max_faces: int = MAX_VOLUME_FACES,
    max_connectivity: int = MAX_VOLUME_CONNECTIVITY,
) -> VolumeSkewnessReport:
    """Compute a reproducible internal-face skewness indicator and adjacency."""
    xyz = np.asarray(points, dtype=np.float64)
    fv = np.asarray(face_vertices)
    fo = np.asarray(face_offsets)
    own = np.asarray(owner)
    nei = np.asarray(neighbour)
    cf = np.asarray(cell_faces)
    co = np.asarray(cell_offsets)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("skewness threshold must be finite and positive")
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError("points must be finite with shape (n_points, 3)")
    for name, array in (("face_vertices", fv), ("face_offsets", fo),
                        ("owner", own), ("neighbour", nei),
                        ("cell_faces", cf), ("cell_offsets", co)):
        if array.ndim != 1 or array.dtype.kind not in "iu":
            raise ValueError(f"{name} must be a one-dimensional integer array")
    n_faces = int(own.size)
    n_cells = int(co.size - 1)
    if n_cells < 1 or n_cells > max_cells or n_faces > max_faces:
        raise ValueError("volume mesh exceeds configured entity limits")
    if fv.size + cf.size > max_connectivity:
        raise ValueError("volume connectivity exceeds configured entry limit")
    if fo.shape != (n_faces + 1,) or nei.shape != (n_faces,):
        raise ValueError("face arrays have inconsistent lengths")
    if co.shape != (n_cells + 1,) or int(fo[0]) != 0 or int(co[0]) != 0:
        raise ValueError("CSR offsets have invalid dimensions or origin")
    if (np.any(fo[1:] < fo[:-1]) or int(fo[-1]) != fv.size
            or np.any(co[1:] < co[:-1]) or int(co[-1]) != cf.size):
        raise ValueError("CSR offsets must be monotonic and terminate at connectivity size")
    if np.any(np.diff(fo) < 3) or np.any(np.diff(co) < 4):
        raise ValueError("each face needs 3 vertices and each cell at least 4 faces")
    if fv.size and (int(fv.min()) < 0 or int(fv.max()) >= xyz.shape[0]):
        raise ValueError("face vertex index is out of range")
    if cf.size and (int(cf.min()) < 0 or int(cf.max()) >= n_faces):
        raise ValueError("cell face index is out of range")
    if np.any(own < 0) or np.any(own >= n_cells):
        raise ValueError("owner cell index is out of range")
    if np.any(nei < -1) or np.any(nei >= n_cells):
        raise ValueError("neighbour cell index is out of range")
    internal = np.flatnonzero(nei >= 0)
    if np.any(own[internal] == nei[internal]):
        raise ValueError("internal face cannot have identical owner and neighbour")
    # Vertex-average cell centroids are robust for arbitrary polyhedra and do
    # not assume convexity or a particular element family.
    cell_centres = np.empty((n_cells, 3), dtype=np.float64)
    for cell_id in range(n_cells):
        face_ids = cf[int(co[cell_id]):int(co[cell_id + 1])]
        point_ids = np.unique(np.concatenate([
            fv[int(fo[fid]):int(fo[fid + 1])] for fid in face_ids
        ]))
        if point_ids.size < 4:
            raise ValueError(f"cell {cell_id} has fewer than four unique vertices")
        cell_centres[cell_id] = xyz[point_ids].mean(axis=0)
    face_scores = np.full(n_faces, np.nan, dtype=np.float64)
    cell_scores = np.zeros(n_cells, dtype=np.float64)
    adjacency = [set() for _ in range(n_cells)]
    for face_id in internal:
        owner_id, neighbour_id = int(own[face_id]), int(nei[face_id])
        adjacency[owner_id].add(neighbour_id)
        adjacency[neighbour_id].add(owner_id)
        polygon = xyz[fv[int(fo[face_id]):int(fo[face_id + 1])]]
        face_centre = polygon.mean(axis=0)
        area_vector = np.sum(np.cross(polygon, np.roll(polygon, -1, axis=0)), axis=0) * 0.5
        area = float(np.linalg.norm(area_vector))
        if not np.isfinite(area) or area <= np.finfo(float).tiny:
            raise ValueError(f"internal face {face_id} is degenerate")
        line = cell_centres[neighbour_id] - cell_centres[owner_id]
        line_norm2 = float(line @ line)
        if line_norm2 <= np.finfo(float).tiny:
            raise ValueError(f"cells adjacent to face {face_id} have coincident centroids")
        projection = float((face_centre - cell_centres[owner_id]) @ line) / line_norm2
        closest = cell_centres[owner_id] + projection * line
        # Dimensionless displacement relative to equivalent face diameter.
        scale = np.sqrt(area)
        score = float(np.linalg.norm(face_centre - closest) / scale)
        face_scores[face_id] = score
        cell_scores[owner_id] = max(cell_scores[owner_id], score)
        cell_scores[neighbour_id] = max(cell_scores[neighbour_id], score)
    if internal.size == 0:
        raise ValueError("mesh has no internal faces; skewness cannot be assessed")
    bad = np.flatnonzero(cell_scores > threshold)
    order = np.argsort(cell_scores)[::-1]
    return VolumeSkewnessReport(
        cell_skewness=cell_scores,
        face_skewness=face_scores,
        cell_neighbours=tuple(tuple(sorted(values)) for values in adjacency),
        worst_cell_ids=tuple(int(i) for i in order[:20]),
        bad_cell_ids=tuple(int(i) for i in bad),
        threshold=float(threshold),
        internal_face_count=int(internal.size),
        max_skewness=float(np.max(face_scores[internal])),
        p95_skewness=float(np.percentile(face_scores[internal], 95)),
    )


def analyze_cfdx_volume_skewness(
    path: str | Path, *, expected_signature: tuple[int, int, int, int] | None = None,
    threshold: float = 0.5,
) -> VolumeSkewnessReport:
    """Load and analyze root CFDX volume topology under strict HDF5 budgets."""
    import h5py

    source = Path(path)
    stat = source.stat()
    signature = (int(stat.st_dev), int(stat.st_ino), int(stat.st_size), int(stat.st_mtime_ns))
    if expected_signature is not None and signature != expected_signature:
        raise RuntimeError("mesh file changed since metadata inspection")
    names = ("points", "face_vertices", "face_offsets", "owner", "neighbour",
             "cell_faces", "cell_offsets")
    with h5py.File(source, "r") as h5:
        arrays = {}
        total_bytes = 0
        for name in names:
            link = h5.get(name, getlink=True)
            dataset = h5.get(name)
            if not isinstance(link, h5py.HardLink) or not isinstance(dataset, h5py.Dataset):
                raise TypeError(f"volume quality requires internal hard-linked {name}")
            if dataset.is_virtual or dataset.external:
                raise ValueError(f"volume quality refuses external dataset {name}")
            if dataset.dtype.kind not in "fiu" or dataset.dtype.kind == "O":
                raise ValueError(f"volume quality requires numeric dataset {name}")
            if dataset.chunks is not None and prod(dataset.chunks) * dataset.dtype.itemsize > 1024 * 1024:
                raise ValueError(f"{name} chunk exceeds the 1 MiB read budget")
            total_bytes += int(dataset.size * dataset.dtype.itemsize)
            arrays[name] = dataset
        if total_bytes > MAX_VOLUME_SOURCE_BYTES:
            raise ValueError("volume topology exceeds the 64 MiB eager-read budget")
        if arrays["points"].ndim != 2 or arrays["points"].shape[1] != 3:
            raise ValueError("points must have shape (n_points, 3)")
        if arrays["owner"].ndim != 1:
            raise ValueError("owner must be one-dimensional")
        if arrays["cell_offsets"].ndim != 1:
            raise ValueError("cell_offsets must be one-dimensional")
        data = {name: dataset[...] for name, dataset in arrays.items()}
    stat_after = source.stat()
    after = (int(stat_after.st_dev), int(stat_after.st_ino), int(stat_after.st_size), int(stat_after.st_mtime_ns))
    if after != signature:
        raise RuntimeError("mesh file changed during volume quality analysis")
    return analyze_volume_skewness(
        data["points"], data["face_vertices"], data["face_offsets"],
        data["owner"], data["neighbour"], data["cell_faces"], data["cell_offsets"],
        threshold=threshold,
    )
