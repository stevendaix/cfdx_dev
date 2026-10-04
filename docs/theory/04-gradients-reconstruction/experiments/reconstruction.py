"""Reference Green--Gauss, LS and WLS gradient reconstructions."""

from __future__ import annotations

import numpy as np

from mesh import Mesh2D


def green_gauss(mesh: Mesh2D, values: np.ndarray) -> np.ndarray:
    """Cell Green--Gauss using arithmetic face interpolation.

    Boundary cells are still computed, but convergence studies should normally
    restrict the error norm to cells with a complete neighbour stencil.
    """

    values = np.asarray(values, dtype=float)
    gradients = np.zeros((len(mesh.centers), 2))
    for p, neighbours in enumerate(mesh.neighbours):
        flux = np.zeros(2)
        for face, nb in enumerate(neighbours):
            if nb >= 0:
                face_value = 0.5 * (values[p] + values[nb])
            else:
                # Reference experiment only: exact face value is not used by
                # the formal interior-cell convergence measurement.
                continue
            flux += face_value * mesh.face_area_vectors[p, face]
        gradients[p] = flux / mesh.volumes[p]
    return gradients


def least_squares(
    mesh: Mesh2D,
    values: np.ndarray,
    weighted: bool = False,
    power: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    gradients = np.full((len(mesh.centers), 2), np.nan)
    condition_numbers = np.full(len(mesh.centers), np.inf)

    for p, neighbours in enumerate(mesh.neighbours):
        valid = neighbours[neighbours >= 0]
        if len(valid) < 2:
            continue

        delta = mesh.centers[valid] - mesh.centers[p]
        rhs = values[valid] - values[p]
        if weighted:
            distance = np.linalg.norm(delta, axis=1)
            weights = 1.0 / np.maximum(distance, np.finfo(float).eps) ** power
        else:
            weights = np.ones(len(valid))

        aw = delta * np.sqrt(weights)[:, None]
        bw = rhs * np.sqrt(weights)
        normal = aw.T @ aw
        condition_numbers[p] = np.linalg.cond(normal)
        if np.linalg.matrix_rank(normal) < 2:
            continue
        gradients[p] = np.linalg.solve(normal, aw.T @ bw)

    return gradients, condition_numbers
