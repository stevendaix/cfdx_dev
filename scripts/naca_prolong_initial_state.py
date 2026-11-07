"""Opt-in coarse-to-fine NACA initial guess; NEVER a qualification result.

Uses KD-tree interpolation to make a strict cell-aligned CFDX CSV from a
previously CONVERGED coarse solution and fine cell centers. Refinement
may accelerate initialization but all nonlinear and physical gates still apply.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
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


def prolong_nested_o_grid(
    coarse: np.ndarray,
    fine: np.ndarray,
    *,
    coarse_surface_n: int,
    fine_surface_n: int,
) -> np.ndarray:
    """Prolong a ring-major CFDX NACA O-grid without crossing wall branches.

    This opt-in interpolation uses four coarse *topological* neighbours and
    a physical-coordinate affine least-squares reconstruction. A KD-tree may
    select cells on the opposite side of the thin trailing-edge wall. Only an
    initial guess is generated: no solution or residual is certified.

    The CFDX mesh generator has (2*n_surface+1) azimuthal cell columns;
    a doubled surface resolution gives m_fine=2*m_coarse-1 because the
    finite-thickness trailing-edge closing segment is not subdivided.
    Radial subdivision must be exactly two. The last azimuthal cell closes
    the seam and remains a single fine cell.
    """
    if (coarse.ndim != 2 or coarse.shape[1] != 6 or
            fine.ndim != 2 or fine.shape[1] != 3 or
            not np.all(np.isfinite(coarse)) or not np.all(np.isfinite(fine))):
        raise ValueError("nested NACA prolongation: invalid or non-finite input")
    if (coarse_surface_n < 2 or fine_surface_n != 2 * coarse_surface_n):
        raise ValueError("nested NACA prolongation: surface count must double")
    mc = 2 * coarse_surface_n + 1
    mf = 2 * fine_surface_n + 1
    if (len(coarse) % mc or len(fine) % mf):
        raise ValueError("nested NACA prolongation: cells are not ring-major")
    nr_c, nr_f = len(coarse) // mc, len(fine) // mf
    if nr_c < 2 or nr_f != 2 * nr_c:
        raise ValueError("nested NACA prolongation: radial count must double")
    # CFDX extrudes one 2D cell layer with identical spanwise centres.
    if (np.ptp(coarse[:, 2]) > 1e-10 or np.ptp(fine[:, 2]) > 1e-10 or
            abs(coarse[0, 2] - fine[0, 2]) > 1e-10):
        raise ValueError("nested NACA prolongation: only planar extrusions supported")

    cr = coarse.reshape(nr_c, mc, 6)
    out = np.empty((len(fine), 6), dtype=np.float64)
    out[:, :3] = fine
    # Conditioning is evaluated on locally dimensionless coordinates, so
    # an extremely small first wall layer is not mistaken for a singular
    # matrix merely because the mesh has a large physical aspect ratio.
    for j in range(nr_f):
        radial = 0.5 * (j + 0.5) - 0.5
        r0 = min(max(int(math.floor(radial)), 0), nr_c - 2)
        r1 = r0 + 1
        for i in range(mf):
            # The last face across the finite TE gap is not subdivided.
            azimuth = 0.5 * i - 0.25 if i < mf - 1 else float(mc - 1)
            a0 = math.floor(azimuth) % mc
            a1 = (a0 + 1) % mc
            sample = np.stack((cr[r0, a0], cr[r0, a1],
                               cr[r1, a0], cr[r1, a1]))
            target = fine[j * mf + i, :2]
            displacement = sample[:, :2] - target
            span = np.max(np.linalg.norm(sample[:, :2] - sample[0, :2], axis=1))
            if not (math.isfinite(span) and span > 0.0):
                raise ValueError("nested NACA prolongation: degenerate coarse stencil")
            if np.min(np.linalg.norm(displacement, axis=1)) > 2.5 * span:
                raise ValueError("nested NACA prolongation: unrelated mesh geometry")
            # Solve q(x,y)=q_f+g_x*(x-x_f)+g_y*(y-y_f) for all
            # momentum/pressure channels. At affine fields this recovers
            # q_f exactly; unstructured inverse-distance averaging does not.
            matrix = np.column_stack((np.ones(4), displacement / span))
            coefficients, _, rank, singular = np.linalg.lstsq(
                matrix, sample[:, 3:], rcond=1e-11)
            if rank != 3 or singular[-1] <= 1e-9 * singular[0]:
                raise ValueError("nested NACA prolongation: rank-deficient stencil")
            fields = coefficients[0]
            if not np.all(np.isfinite(fields)):
                raise ValueError("nested NACA prolongation: non-finite initial guess")
            # A pathological extrapolation should fail instead of silently
            # creating a massive initial pressure or velocity spike.
            channel_scale = np.maximum(1.0, np.max(np.abs(sample[:, 3:]), axis=0))
            if np.any(np.abs(fields) > 10.0 * channel_scale):
                raise ValueError("nested NACA prolongation: unstable extrapolation")
            out[j * mf + i, 3:] = fields
    return out


def create(
    coarse_path: Path, fine_path: Path, output: Path,
    *, method: str = "nearest-idw", coarse_surface_n: int | None = None,
    fine_surface_n: int | None = None,
) -> int:
    coarse = read_csv(coarse_path, ("cell", "x", "y", "z", "ux", "uy", "p"))
    fine = read_csv(fine_path, ("cell", "x", "y", "z"))
    if method == "nested-wls":
        if coarse_surface_n is None or fine_surface_n is None:
            raise ValueError("nested NACA prolongation needs both surface counts")
        data = prolong_nested_o_grid(
            coarse, fine, coarse_surface_n=coarse_surface_n,
            fine_surface_n=fine_surface_n)
    elif method == "nearest-idw":
        if coarse_surface_n is not None or fine_surface_n is not None:
            raise ValueError("surface counts require nested-wls")
        data = prolong(coarse, fine)
    else:
        raise ValueError("unknown NACA prolongation method")
    # The output is an initial guess, never a solution certificate. Do not
    # destroy the source checkpoint or the target mesh centres when someone
    # accidentally reuses an input path. Publish only a complete CSV.
    if output.resolve() in (coarse_path.resolve(), fine_path.resolve()):
        raise ValueError("NACA prolongation output must differ from input files")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", newline="", encoding="utf-8", suffix=".incomplete",
            prefix=output.name + ".", dir=output.parent, delete=False
        ) as stream:
            temporary = Path(stream.name)
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(("cell", "x", "y", "z", "ux", "uy", "p"))
            for i, row in enumerate(data):
                writer.writerow((i, *[format(float(x), ".17g") for x in row]))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, output)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    print(f"NACA_PROLONG: INITIAL_GUESS_ONLY method={method} "
          f"coarse={len(coarse)} fine={len(fine)}")
    return len(fine)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine-centers", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--method", choices=("nearest-idw", "nested-wls"),
                        default="nearest-idw")
    parser.add_argument("--coarse-surface-n", type=int)
    parser.add_argument("--fine-surface-n", type=int)
    args = parser.parse_args()
    try:
        create(args.coarse, args.fine_centers, args.output,
               method=args.method, coarse_surface_n=args.coarse_surface_n,
               fine_surface_n=args.fine_surface_n)
    except (ValueError, OSError, OverflowError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
