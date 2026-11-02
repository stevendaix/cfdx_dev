"""Opt-in coarse-to-fine NACA initial guess; NEVER a qualification result.

Uses KD-tree interpolation to make a strict cell-aligned CFDX CSV from a
previously CONVERGED coarse solution and fine cell centers. Refinement
may accelerate initialization but all nonlinear and physical gates still apply.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def read_csv(path: Path, required: tuple[str, ...]) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != required:
            raise ValueError(f"bad CSV columns in {path}")
        rows = list(reader)
    if not rows or any(int(row["cell"]) != i for i, row in enumerate(rows)):
        raise ValueError("invalid or unordered cell ids")
    data = np.asarray([[float(row[key]) for key in required[1:]] for row in rows])
    if not np.all(np.isfinite(data)):
        raise ValueError("non-finite state or geometry")
    return data


def prolong(coarse: np.ndarray, fine: np.ndarray, neighbors: int = 6) -> np.ndarray:
    if coarse.ndim != 2 or coarse.shape[1] != 6:
        raise ValueError("coarse field must be an N x 6 array")
    if fine.ndim != 2 or fine.shape[1] != 3:
        raise ValueError("fine mesh must be an M x 3 array")
    if not (4 <= neighbors <= 16) or coarse.shape[0] < neighbors:
        raise ValueError("insufficient neighbors for prolongation")
    if not np.all(np.isfinite(coarse)) or not np.all(np.isfinite(fine)):
        raise ValueError("non-finite input")
    # The NACA case is an extruded 2-D flow; use physical x/y only.
    tree = cKDTree(coarse[:, :2])
    distance, index = tree.query(fine[:, :2], k=neighbors)
    nearest = distance[:, 0]
    # Do not extrapolate to a different computational domain.
    scale = np.sqrt(np.ptp(coarse[:, 0]) ** 2 + np.ptp(coarse[:, 1]) ** 2)
    if not math.isfinite(scale) or not scale > 0:
        raise ValueError("degenerate coarse coordinates")
    envelope = 0.02 * scale
    if (np.any(fine[:, :2] < coarse[:, :2].min(axis=0) - envelope)
            or np.any(fine[:, :2] > coarse[:, :2].max(axis=0) + envelope)):
        raise ValueError("fine centers outside the coarse domain envelope")
    near = nearest < 1e-14 * max(scale, 1)
    weight = 1 / np.maximum(distance, 1e-14 * max(scale, 1)) ** 2
    weight /= weight.sum(axis=1)[:, None]
    fields = np.einsum("ij,ijk->ik", weight, coarse[index, 3:])
    if np.any(near):
        fields[near] = coarse[index[near, 0], 3:]
    if not np.all(np.isfinite(fields)):
        raise ValueError("non-finite interpolated field")
    return np.column_stack((fine, fields))


def create(coarse_path: Path, fine_path: Path, output: Path) -> int:
    coarse = read_csv(coarse_path, ("cell", "x", "y", "z", "ux", "uy", "p"))
    fine = read_csv(fine_path, ("cell", "x", "y", "z"))
    data = prolong(coarse, fine)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("cell", "x", "y", "z", "ux", "uy", "p"))
        for i, row in enumerate(data):
            writer.writerow((i, *[format(float(x), ".17g") for x in row]))
    print(f"NACA_PROLONG: INITIAL_GUESS_ONLY coarse={len(coarse)} fine={len(fine)}")
    return len(fine)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine-centers", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        create(args.coarse, args.fine_centers, args.output)
    except (ValueError, OSError, OverflowError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
