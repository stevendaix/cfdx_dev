#!/usr/bin/env python3
"""Import the single-block NASA/TMR NACA0012 PLOT3D C-grid.

The j=0 boundary contains the airfoil AND two coincident wake cuts.  The
wake cuts are joined as internal faces; only the closed airfoil is a wall.
This is an input/topology conversion, not a physical qualification.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path
import re

import h5py
import numpy as np

from generate_2d_validation_meshes import write_mesh

TMR_URL = "https://tmbwg.github.io/turbmodels/naca0012_grids.html"
TMR_ARCHIVE_URL = "https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip"


def read_plot3d_2d(path: Path):
    """Read NASA's formatted multigrid PLOT3D, with one 2-D block.

    The first record is nbl, not ni.  Compressed files in the advertised
    grid family use the same format and are read without manual unpacking.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        tokens = stream.read().replace(",", " ").split()
    if len(tokens) < 3:
        raise ValueError(f"PLOT3D block header is missing: {path}")
    nblocks = int(tokens[0])
    if nblocks != 1:
        raise ValueError(f"only single-block NASA/TMR PLOT3D is supported, got {nblocks}")
    ni, nj = int(tokens[1]), int(tokens[2])
    if ni < 5 or nj < 2:
        raise ValueError(f"invalid NASA/TMR C-grid dimensions: {ni}x{nj}")
    expected = 3 + 2 * ni * nj
    if len(tokens) != expected:
        raise ValueError(f"expected {expected} numeric tokens for {ni}x{nj}, got {len(tokens)}")
    values = np.asarray(
        [float(value.replace("D", "E").replace("d", "e")) for value in tokens[3:]],
        dtype=np.float64,
    )
    if not np.isfinite(values).all():
        raise ValueError("NASA/TMR grid coordinates must be finite")
    return values[:ni * nj].reshape((nj, ni)), values[ni * nj:].reshape((nj, ni))


def convert(x, y):
    """Join the coincident j=0 wake and extrude with EMPTY front/back.

    NASA's neutral map declares a reversed one-to-one connection between
    the leading and trailing index intervals on j=0.  Determine that
    interval from coincident endpoint pairs, including the trailing-edge
    pair.  No wall-index ratio is guessed from the filename.  Reject a
    different topology instead of manufacturing a wall in its wake.
    """
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if x.ndim != 2 or x.shape != y.shape or min(x.shape) < 2:
        raise ValueError("C-grid x/y arrays must have the same two-dimensional shape")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("C-grid coordinates must be finite")
    nj, ni = x.shape
    # Absolute tolerance is in the unit-chord NASA grid coordinates.  It
    # is deliberately independent of the distant farfield radius.
    paired_points = 0
    for i in range(ni // 2):
        opposite = ni - 1 - i
        if max(abs(x[0, i] - x[0, opposite]),
               abs(y[0, i] - y[0, opposite])) > 1e-12:
            break
        paired_points += 1
    first_wall = paired_points - 1
    last_wall = ni - paired_points
    if paired_points < 2 or last_wall - first_wall < 3:
        raise ValueError("expected a reversed coincident j=0 wake and a nondegenerate airfoil")

    # The supplied NASA family is normalized to unit chord.  Check the
    # independently known geometry before accepting the detected seam: a
    # partly corrupted wake must not be reclassified as extra solid wall.
    wall_x = x[0, first_wall:last_wall + 1]
    wall_y = y[0, first_wall:last_wall + 1]
    geometry_tolerance = 1e-6  # admits the documented original TMR TE typo.
    if (abs(wall_x[0] - 1.0) > geometry_tolerance or
            abs(wall_y[0]) > geometry_tolerance or
            abs(wall_x.min()) > geometry_tolerance or
            wall_x.max() > 1.0 + geometry_tolerance or
            (np.diff(x[0, :paired_points]) >= 0.0).any()):
        raise ValueError("matched wake must terminate at the unit-chord NACA trailing edge")

    # Merge only the explicitly matched wake vertices.  This gives each
    # wake face two owners through write_mesh's ordinary connectivity
    # matching; it is not a periodic/no-slip physical boundary.
    representatives = np.arange(ni * nj)
    for i in range(paired_points):
        representatives[ni - 1 - i] = i
    unique, inverse = np.unique(representatives, return_inverse=True)
    nplane = len(unique)
    flat_x, flat_y = x.ravel(), y.ravel()
    points = [(float(flat_x[i]), float(flat_y[i]), z)
              for z in (0.0, 1.0) for i in unique]

    def p(j, i, layer=0):
        return int(inverse[j * ni + i]) + layer * nplane

    cells = []
    for j in range(nj - 1):
        for i in range(ni - 1):
            lower = (p(j, i), p(j, i + 1), p(j + 1, i + 1), p(j + 1, i))
            cell = lower + tuple(v + nplane for v in lower)
            if len(set(cell)) != 8:
                raise ValueError("wake stitching created a degenerate hexahedron")
            cells.append(cell)
    airfoil = [(p(0, i + 1), p(0, i), p(0, i, 1), p(0, i + 1, 1))
               for i in range(first_wall, last_wall)]
    farfield = [(p(nj - 1, i), p(nj - 1, i + 1),
                 p(nj - 1, i + 1, 1), p(nj - 1, i, 1))
                for i in range(ni - 1)]
    # i=0 and i=ni-1 are downstream farfield boundaries in the neutral
    # map.  The actual wake connection is along j=0, handled above.
    farfield += [(p(j, 0), p(j + 1, 0), p(j + 1, 0, 1), p(j, 0, 1))
                 for j in range(nj - 1)]
    farfield += [(p(j + 1, ni - 1), p(j, ni - 1),
                  p(j, ni - 1, 1), p(j + 1, ni - 1, 1))
                 for j in range(nj - 1)]
    front = [tuple(cell[:4]) for cell in cells]
    back = [tuple(reversed(cell[4:])) for cell in cells]
    return points, cells, {"airfoil": airfoil, "farfield": farfield,
                           "front": front, "back": back}


def infer_dimensions(path: Path):
    match = re.search(r"_(\d+)-(\d+)\.p2dfmt(?:\.gz)?$", path.name)
    if not match:
        raise ValueError(f"cannot infer TMR grid dimensions from {path.name}")
    return int(match.group(1)), int(match.group(2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ni_name, nj_name = infer_dimensions(args.input)
    x, y = read_plot3d_2d(args.input)
    nj, ni = x.shape
    if (ni, nj) != (ni_name, nj_name):
        raise ValueError(f"filename declares {ni_name}x{nj_name}, file contains {ni}x{nj}")
    points, cells, boundaries = convert(x, y)
    write_mesh(args.output, points, cells, boundaries)
    with h5py.File(args.output, "r+") as h5:
        h5.attrs["validation_source"] = "NASA TMR NACA0012"
        h5.attrs["validation_source_url"] = TMR_URL
        h5.attrs["validation_archive_url"] = TMR_ARCHIVE_URL
        h5.attrs["validation_source_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
        h5.attrs["validation_grid"] = f"{ni}x{nj}"
        h5.attrs["validation_topology"] = "C-grid; j=0 reversed wake stitched"
        h5.attrs["validation_wall_faces"] = len(boundaries["airfoil"])
        h5.attrs["validation_wake_pairs"] = (ni - 1 - len(boundaries["airfoil"])) // 2
        h5.attrs["validation_physics"] = "N9-S8 laminar incompressible campaign"
        h5.attrs["validation_status"] = "imported-provenance-only"
    print(f"Imported {args.input} -> {args.output} ({ni}x{nj}); "
          f"airfoil_faces={len(boundaries['airfoil'])}; physical qualification NOT RUN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
