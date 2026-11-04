"""Bounded geometric quality diagnostics for explicitly selected mesh surfaces.

These diagnostics describe polygon geometry only. They are not substitutes for
volume-cell quality metrics such as non-orthogonality, skewness, Jacobian or
cell aspect ratio, which require complete cell geometry and adjacency.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import prod
from pathlib import Path
from typing import Any

import numpy as np


def _file_signature(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return (
        int(stat.st_dev),
        int(stat.st_ino),
        int(stat.st_size),
        int(stat.st_mtime_ns),
    )


def analyze_cfdx_surface_file(
    path: str | Path,
    *,
    expected_signature: tuple[int, int, int, int] | None = None,
    max_source_bytes: int = 32 * 1024 * 1024,
    max_block_bytes: int = 1024 * 1024,
) -> SurfaceQualityReport:
    """Read and analyze the supported root surface topology under strict budgets."""
    import h5py

    if max_source_bytes <= 0 or max_block_bytes <= 0:
        raise ValueError("surface quality byte budgets must be positive")
    source = Path(path)
    signature_before = _file_signature(source)
    if expected_signature is not None and signature_before != expected_signature:
        raise RuntimeError("mesh file changed since metadata inspection")
    with h5py.File(source, "r") as h5:
        datasets = {}
        for name in ("points", "face_vertices", "face_offsets"):
            link = h5.get(name, getlink=True)
            dataset = h5.get(name)
            if not isinstance(link, h5py.HardLink) or not isinstance(
                dataset, h5py.Dataset
            ):
                raise ValueError(
                    f"surface quality requires an internal hard-linked {name} dataset"
                )
            if dataset.is_virtual or dataset.external:
                raise ValueError(f"surface quality does not read external {name}")
            if dataset.chunks is not None:
                chunk_bytes = prod(dataset.chunks) * dataset.dtype.itemsize
                if chunk_bytes > max_block_bytes:
                    raise ValueError(
                        f"{name} chunk exceeds the {max_block_bytes:,}-byte block budget"
                    )
            datasets[name] = dataset

        points, vertices, offsets = (
            datasets["points"],
            datasets["face_vertices"],
            datasets["face_offsets"],
        )
        if points.ndim != 2 or points.shape[1] != 3 or points.dtype.kind != "f":
            raise ValueError("points must be a floating-point (n_points, 3) array")
        if vertices.ndim != 1 or vertices.dtype.kind not in "iu":
            raise ValueError("face_vertices must be a one-dimensional integer array")
        if offsets.ndim != 1 or offsets.dtype.kind not in "iu":
            raise ValueError("face_offsets must be a one-dimensional integer array")
        face_count = int(offsets.shape[0] - 1)
        if face_count < 0 or face_count > 100_000:
            raise ValueError("surface face count exceeds the 100k quality limit")
        if vertices.shape[0] > 2_000_000:
            raise ValueError("surface connectivity exceeds the 2M-entry quality limit")
        source_bytes = sum(
            prod(dataset.shape) * dataset.dtype.itemsize
            for dataset in datasets.values()
        )
        if source_bytes > max_source_bytes:
            raise ValueError(
                f"surface topology requires {source_bytes:,} bytes; "
                f"budget is {max_source_bytes:,}"
            )
        point_values = points[...]
        face_values = vertices[...]
        face_offsets = offsets[...]
    if _file_signature(source) != signature_before:
        raise RuntimeError("mesh file changed during surface quality analysis")
    return analyze_surface_quality(
        point_values,
        face_values,
        face_offsets,
    )


@dataclass(frozen=True)
class SurfaceQualityReport:
    """Summary of geometric quality indicators for polygonal surface faces."""

    face_count: int
    valid_face_count: int
    degenerate_face_ids: tuple[int, ...]
    zero_length_edge_face_ids: tuple[int, ...]
    face_area_min: float | None
    face_area_max: float | None
    edge_length_ratio_min: float | None
    edge_length_ratio_max: float | None
    edge_length_ratio_p95: float | None
    warpage_ratio_max: float | None
    warpage_ratio_p95: float | None
    worst_warpage_face_ids: tuple[int, ...]


def analyze_surface_quality(
    points: Any,
    face_vertices: Any,
    face_offsets: Any,
    *,
    max_faces: int = 100_000,
    max_connectivity_entries: int = 2_000_000,
    max_reported_faces: int = 20,
) -> SurfaceQualityReport:
    """Compute bounded face area, edge-ratio and plane-warpage diagnostics.

    The edge-length ratio is max(edge length) / min(edge length), and is
    intentionally not called a cell aspect ratio. Warpage is the maximum
    absolute distance from face vertices to their least-squares plane,
    divided by the face's longest edge. Degeneracy uses a scale-relative
    floating-point threshold rather than a mesh-unit-specific tolerance.

    Arrays must already be in memory; callers own file-I/O byte budgets.
    """
    if max_faces <= 0 or max_connectivity_entries <= 0 or max_reported_faces < 0:
        raise ValueError("surface quality limits must be positive")
    coordinates = np.asarray(points)
    vertices = np.asarray(face_vertices)
    offsets = np.asarray(face_offsets)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("points must have shape (n_points, 3)")
    if coordinates.dtype.kind not in "fiu":
        raise ValueError("point coordinates must be real numeric values")
    if not np.isfinite(coordinates).all():
        raise ValueError("point coordinates must all be finite")
    if vertices.ndim != 1 or vertices.dtype.kind not in "iu":
        raise ValueError("face_vertices must be a one-dimensional integer array")
    if offsets.ndim != 1 or offsets.dtype.kind not in "iu":
        raise ValueError("face_offsets must be a one-dimensional integer array")
    face_count = int(offsets.size - 1)
    if face_count < 0 or face_count > max_faces:
        raise ValueError("surface face count exceeds the configured quality limit")
    if vertices.size > max_connectivity_entries:
        raise ValueError("surface connectivity exceeds the configured quality limit")
    if offsets.size == 0 or int(offsets[0]) != 0:
        raise ValueError("face_offsets must start at zero")
    if np.any(offsets[1:] < offsets[:-1]) or int(offsets[-1]) != vertices.size:
        raise ValueError(
            "face_offsets must be monotonic and terminate at connectivity size"
        )
    if vertices.size and (
        int(vertices.min()) < 0 or int(vertices.max()) >= coordinates.shape[0]
    ):
        raise ValueError("face_vertices contains an out-of-range point index")
    if np.any(np.diff(offsets) < 3):
        raise ValueError("every surface face must contain at least three vertices")

    areas: list[float] = []
    edge_ratios: list[float] = []
    warpages: list[float] = []
    worst_ids: list[tuple[float, int]] = []
    degenerate: list[int] = []
    zero_edges: list[int] = []

    for face_id in range(face_count):
        start, end = int(offsets[face_id]), int(offsets[face_id + 1])
        polygon = coordinates[vertices[start:end]].astype(np.float64, copy=False)
        edges = np.roll(polygon, -1, axis=0) - polygon
        edge_lengths = np.linalg.norm(edges, axis=1)
        longest = float(np.max(edge_lengths))
        shortest = float(np.min(edge_lengths))
        if longest == 0.0 or shortest <= np.finfo(np.float64).eps * longest:
            zero_edges.append(face_id)

        centered = polygon - np.mean(polygon, axis=0)
        area_vector = np.sum(
            np.cross(centered, np.roll(centered, -1, axis=0)), axis=0
        )
        area = 0.5 * float(np.linalg.norm(area_vector))
        area_tolerance = 64.0 * np.finfo(np.float64).eps * longest * longest
        if longest == 0.0 or area <= area_tolerance:
            degenerate.append(face_id)
            continue

        areas.append(area)
        edge_ratios.append(longest / shortest)
        try:
            _, _, right_vectors = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f"could not determine face plane for face {face_id}") from exc
        warpage = float(
            np.max(np.abs(centered @ right_vectors[-1])) / longest
        )
        warpages.append(warpage)
        worst_ids.append((warpage, face_id))

    if not areas:
        return SurfaceQualityReport(
            face_count=face_count,
            valid_face_count=0,
            degenerate_face_ids=tuple(degenerate),
            zero_length_edge_face_ids=tuple(zero_edges),
            face_area_min=None,
            face_area_max=None,
            edge_length_ratio_min=None,
            edge_length_ratio_max=None,
            edge_length_ratio_p95=None,
            warpage_ratio_max=None,
            warpage_ratio_p95=None,
            worst_warpage_face_ids=(),
        )

    worst_ids.sort(reverse=True)
    return SurfaceQualityReport(
        face_count=face_count,
        valid_face_count=len(areas),
        degenerate_face_ids=tuple(degenerate),
        zero_length_edge_face_ids=tuple(zero_edges),
        face_area_min=float(np.min(areas)),
        face_area_max=float(np.max(areas)),
        edge_length_ratio_min=float(np.min(edge_ratios)),
        edge_length_ratio_max=float(np.max(edge_ratios)),
        edge_length_ratio_p95=float(np.percentile(edge_ratios, 95)),
        warpage_ratio_max=float(np.max(warpages)),
        warpage_ratio_p95=float(np.percentile(warpages, 95)),
        worst_warpage_face_ids=tuple(
            face_id for _, face_id in worst_ids[:max_reported_faces]
        ),
    )
