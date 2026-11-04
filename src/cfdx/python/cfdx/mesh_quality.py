"""Bounded geometric quality diagnostics for explicitly selected mesh surfaces.

These diagnostics describe polygon geometry only. They are not substitutes for
volume-cell quality metrics such as non-orthogonality, skewness, Jacobian or
cell aspect ratio, which require complete cell geometry and adjacency.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


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

        area_vector = np.sum(
            np.cross(polygon, np.roll(polygon, -1, axis=0)), axis=0
        )
        area = 0.5 * float(np.linalg.norm(area_vector))
        area_tolerance = 64.0 * np.finfo(np.float64).eps * longest * longest
        if longest == 0.0 or area <= area_tolerance:
            degenerate.append(face_id)
            continue

        areas.append(area)
        edge_ratios.append(longest / shortest)
        centered = polygon - np.mean(polygon, axis=0)
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
