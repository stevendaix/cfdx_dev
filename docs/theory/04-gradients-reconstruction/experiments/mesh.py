"""Structured skewed quadrilateral mesh used by the N2 theory experiments.

The mesh is deliberately independent of CFDX implementation code. It is a
small analytical reference model for studying reconstruction algorithms.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Mesh2D:
    vertices: np.ndarray
    cells: np.ndarray
    centers: np.ndarray
    volumes: np.ndarray
    neighbours: np.ndarray
    face_centers: np.ndarray
    face_area_vectors: np.ndarray


def make_skewed_mesh(n: int, skew: float = 0.20) -> Mesh2D:
    if n < 4:
        raise ValueError("n must be at least 4")

    xi = np.linspace(0.0, 1.0, n + 1)
    eta = np.linspace(0.0, 1.0, n + 1)
    vertices = np.empty((n + 1, n + 1, 2), dtype=float)
    for j, e in enumerate(eta):
        for i, x in enumerate(xi):
            vertices[j, i] = (
                x + skew * np.sin(np.pi * e) * x * (1.0 - x),
                e + 0.08 * skew * np.sin(np.pi * x) * e * (1.0 - e),
            )

    cell_vertices = []
    centers = []
    volumes = []
    neighbours = []
    face_centers = []
    face_area_vectors = []

    def vid(j: int, i: int) -> int:
        return j * (n + 1) + i

    flat_vertices = vertices.reshape((-1, 2))

    for j in range(n):
        for i in range(n):
            ids = [vid(j, i), vid(j, i + 1), vid(j + 1, i + 1), vid(j + 1, i)]
            pts = flat_vertices[ids]
            cell_vertices.append(ids)
            centers.append(np.mean(pts, axis=0))

            area = 0.5 * abs(
                np.dot(pts[:, 0], np.roll(pts[:, 1], -1))
                - np.dot(pts[:, 1], np.roll(pts[:, 0], -1))
            )
            volumes.append(area)

            nb = []
            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                ii, jj = i + di, j + dj
                nb.append(jj * n + ii if 0 <= ii < n and 0 <= jj < n else -1)
            neighbours.append(nb)

            fc = []
            fa = []
            for a, b, nb_idx in (
                (3, 0, nb[0]),  # left
                (1, 2, nb[1]),  # right
                (0, 1, nb[2]),  # bottom
                (2, 3, nb[3]),  # top
            ):
                p0, p1 = pts[a], pts[b]
                fc.append(0.5 * (p0 + p1))
                edge = p1 - p0
                # Outward normal for CCW cell ordering.
                fa.append(np.array([edge[1], -edge[0]]))
            face_centers.append(fc)
            face_area_vectors.append(fa)

    return Mesh2D(
        vertices=flat_vertices,
        cells=np.asarray(cell_vertices, dtype=int),
        centers=np.asarray(centers),
        volumes=np.asarray(volumes),
        neighbours=np.asarray(neighbours, dtype=int),
        face_centers=np.asarray(face_centers),
        face_area_vectors=np.asarray(face_area_vectors),
    )


def interior_indices(mesh: Mesh2D) -> np.ndarray:
    return np.flatnonzero(np.all(mesh.neighbours >= 0, axis=1))
