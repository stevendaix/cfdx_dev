"""Geometry-only regression for NACA O-grid symmetry and 3-D extrusion.

Not a substitute for Cd/Cp/Cf/RANS/CHT full physical validation.
"""

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "generate_2d_validation_meshes.py"
SPEC = importlib.util.spec_from_file_location(
    "cfdx_generate_2d_validation_meshes", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MESHGEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MESHGEN)


def xy_area(points, face):
    v = np.asarray([points[i][:2] for i in face])
    return 0.5 * float(
        np.dot(v[:, 0], np.roll(v[:, 1], -1)) - np.dot(v[:, 1], np.roll(v[:, 0], -1))
    )


def verify(n, nr, layers):
    points, cells, patches = MESHGEN.naca_o_grid(n, nr, n_spanwise=layers, span=2.0)
    pts = np.asarray(points)
    m = 2 * n + 1
    assert len(cells) == m * nr * layers
    assert pts.shape == ((nr + 1) * m * (layers + 1), 3)
    assert len(patches["front"]) == len(patches["back"]) == m * nr
    assert len(patches["airfoil"]) == len(patches["farfield"]) == m * layers
    for k in range(layers + 1):
        for j in (0, nr):
            ring = pts[(k * (nr + 1) + j) * m : (k * (nr + 1) + j + 1) * m]
            mismatch = float(np.max(np.abs(ring[:, 0] - ring[::-1, 0])))
            mismatch = max(mismatch, float(np.max(np.abs(ring[:, 1] + ring[::-1, 1]))))
            assert mismatch < 5e-12, (n, nr, layers, j, mismatch)

    # Positive cell volumes on the mapped hexahedral grid.
    min_volume = float("inf")
    orientation = None
    for cell in cells:
        a = xy_area(points, cell[:4])
        dz = points[cell[4]][2] - points[cell[0]][2]
        volume = a * dz
        if orientation is None:
            orientation = np.sign(volume)
        assert orientation * volume > 0.0 and np.isfinite(volume), (cell, volume)
        min_volume = min(min_volume, abs(volume))
    assert min_volume > 0.0
    # Lateral face area summed through layers must equal the perimeter x span.
    contour = pts[:m, :2]
    perimeter = float(
        np.linalg.norm(contour - np.roll(contour, -1, axis=0), axis=1).sum()
    )
    area = sum(
        np.linalg.norm(pts[f[1]] - pts[f[0]]) * np.linalg.norm(pts[f[3]] - pts[f[0]])
        for f in patches["airfoil"]
    )
    assert abs(area / (perimeter * 2.0) - 1.0) < 1e-12
    print(
        f"NACA_MESH: PASS n={n} radial={nr} layers={layers} "
        f"cells={len(cells)} min_vol={min_volume:.6e} wall_area={area:.8f}"
    )


if __name__ == "__main__":
    for n, nr in ((16, 8), (64, 24), (128, 48)):
        for nz in (1, 2):
            verify(n, nr, nz)
