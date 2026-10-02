"""Deterministic checks for the standalone N2 reference algorithms."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from linear_field import exact_gradient, field
from mesh import interior_indices, make_skewed_mesh
from reconstruction import green_gauss, least_squares


def main() -> None:
    mesh = make_skewed_mesh(16)
    values = field(mesh.centers[:, 0], mesh.centers[:, 1])
    exact = np.broadcast_to(exact_gradient(), (len(mesh.centers), 2))
    idx = interior_indices(mesh)

    gg = green_gauss(mesh, values)[idx]
    assert np.max(np.abs(gg - exact[idx])) < 1.0e-12

    ls, cond = least_squares(mesh, values)
    wls, wcond = least_squares(mesh, values, weighted=True)
    assert np.max(np.abs(ls[idx] - exact[idx])) < 1.0e-12
    assert np.max(np.abs(wls[idx] - exact[idx])) < 1.0e-12
    assert np.all(np.isfinite(cond[idx]))
    assert np.all(np.isfinite(wcond[idx]))

    print("N2 reference self-check: PASS")


if __name__ == "__main__":
    main()
