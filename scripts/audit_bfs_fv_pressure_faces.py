"""Fail-closed bridge from a captured CFDX pressure CSR to *real mesh faces*.

This audits ONLY the implicit, orthogonal pressure-correction operator. The
matrix coefficients come from CFDX's actual post-BC pressure assembly, while
face owner/neighbour topology comes from the original CFDX HDF5 mesh. A high/low
pressure correction is used without first rounding it to one float64.

It cannot certify the full Rhie--Chow flux (including predictor, nonorthogonal
terms, SIMPLEC response and velocity update), physical convergence, or SST.
The captured CSR has no cryptographic mesh identity: exact graph matching is
verified but source/mesh provenance remains an additional required contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
from pathlib import Path

import h5py
import numpy as np
from audit_bfs_compensated_pressure import compensated_residual, two_sum_update
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu
from seed_bfs_pressure_from_capture import read_capture


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_capture_checksum(
    path: Path, matrix: csr_matrix, rhs: np.ndarray, expected_footer: str
) -> None:
    """Reproduce the native cfdx linear capture's FNV64 bit-pattern checksum.

    A footer copied from another file is not an integrity check. The native
    reader hashes metadata, CSR offsets/indices/values, RHS and initial guess,
    each as an IEEE754 or integer 64-bit little-endian word.
    """
    with path.open(encoding="ascii") as stream:
        if stream.readline().strip() != "CFDX_LINEAR_REPLAY_V1":
            raise ValueError("invalid capture header")
        problem, nullspace, scaling = map(int, stream.readline().split())
        nrows, ncols, nnz, budget = map(int, stream.readline().split())
        tolerance_bits = int(stream.readline().strip(), 16)
    a, b, _ = read_capture(path)
    if a.shape != (nrows, ncols) or a.nnz != nnz or not np.array_equal(b, rhs):
        raise ValueError("capture dimensions/data changed during verification")
    words = [nrows, ncols, nnz, problem, nullspace, scaling, budget, tolerance_bits]
    words.extend(int(value) for value in a.indptr)
    for column, bits in zip(a.indices, a.data.view(np.uint64)):
        words.append(int(column))
        words.append(int(bits))
    # read_capture does not expose x0, so recover its 64-bit representation
    # from the original capture's end-of-file data lines.
    with path.open(encoding="ascii") as stream:
        for _ in range(4 + (nrows + 1) + nnz):
            stream.readline()
        for i in range(nrows):
            row = stream.readline().split()
            if len(row) != 2:
                raise ValueError("truncated captured initial guess")
            rhs_bits = int(row[0], 16)
            if rhs_bits != struct.unpack("<Q", struct.pack("<d", float(rhs[i])))[0]:
                raise ValueError("captured RHS checksum/data mismatch")
            words.extend((rhs_bits, int(row[1], 16)))
    state = 14695981039346656037
    mask = (1 << 64) - 1
    for word in words:
        for shift in range(0, 64, 8):
            state = ((state ^ ((word >> shift) & 255)) * 1099511628211) & mask
    if f"{state:016x}" != expected_footer:
        raise ValueError("captured CSR checksum mismatch")


def read_split(
    path: Path, count: int, expected_checksum: str
) -> tuple[np.ndarray, np.ndarray]:
    with path.open(encoding="ascii") as stream:
        if stream.readline().strip() != "CFDX_PRESSURE_SPLIT_V1":
            raise ValueError("invalid CFDX pressure-split magic")
        if stream.readline().strip() != f"{count} {expected_checksum}":
            raise ValueError("pressure split belongs to a different captured CSR")
        high = np.empty(count, dtype=np.float64)
        low = np.empty(count, dtype=np.float64)
        for row in range(count):
            tokens = stream.readline().split()
            if len(tokens) != 2:
                raise ValueError("truncated or malformed high/low row")
            high[row], low[row] = (float.fromhex(word) for word in tokens)
        if (
            stream.read().strip()
            or not np.isfinite(high).all()
            or not np.isfinite(low).all()
        ):
            raise ValueError("extra or nonfinite pressure-split data")
    return high, low


def read_mesh(
    path: Path, pressure_patch: str
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as mesh:
        required = ("owner", "neighbour", "patch_face_offsets", "patch_face_ids")
        if not all(key in mesh for key in required):
            raise ValueError("mesh does not contain native FV face ownership")
        n = int(mesh.attrs["n_cells"])
        owner = mesh["owner"][:].astype(np.int64)
        neighbour = mesh["neighbour"][:].astype(np.int64)
        offsets = mesh["patch_face_offsets"][:]
        face_ids = mesh["patch_face_ids"][:].astype(np.int64)
        names = str(mesh.attrs["boundary_patches"]).split(";")
        patches = [item.split(":")[0] for item in names]
        if len(patches) + 1 != len(offsets) or patches.count(pressure_patch) != 1:
            raise ValueError("missing or ambiguous pressure Dirichlet outlet patch")
        i = patches.index(pressure_patch)
        outlet_faces = face_ids[int(offsets[i]) : int(offsets[i + 1])]
    nf = len(owner)
    if (
        n < 2
        or len(neighbour) != nf
        or len(outlet_faces) == 0
        or np.any(owner < 0)
        or np.any(owner >= n)
        or np.any(neighbour < -1)
        or np.any(neighbour >= n)
        or np.any((neighbour >= 0) & (neighbour == owner))
        or np.any(outlet_faces < 0)
        or np.any(outlet_faces >= nf)
        or len(np.unique(outlet_faces)) != len(outlet_faces)
        or np.any(neighbour[outlet_faces] != -1)
    ):
        raise ValueError("invalid FV topology or pressure boundary face list")
    return n, owner, neighbour, outlet_faces


def verify_face_operator(
    matrix: csr_matrix,
    rhs: np.ndarray,
    high: np.ndarray,
    low: np.ndarray,
    owner: np.ndarray,
    neighbour: np.ndarray,
    outlet_faces: np.ndarray,
    *,
    target_tolerance: float,
) -> dict:
    """Check CSR/real-mesh edge bijection and reconstruct pressure flux divergence.

    The only supported boundary is one Dirichlet-pressure face per anchored
    cell. No coefficient is guessed or distributed across multiple faces.
    Mixed/nonorthogonal/Schur terms cannot be inferred from CSR and are not
    claimed here.
    """
    a = matrix.tocsr(copy=True)
    a.sum_duplicates()
    a.sort_indices()
    n = a.shape[0]
    if (
        n < 2
        or a.shape != (n, n)
        or len(rhs) != n
        or len(high) != n
        or len(low) != n
        or len(owner) != len(neighbour)
        or not (np.isfinite(target_tolerance) and target_tolerance > 0)
        or not np.isfinite(a.data).all()
        or not np.isfinite(rhs).all()
        or not np.isfinite(high).all()
        or not np.isfinite(low).all()
        or np.any(owner < 0)
        or np.any(owner >= n)
        or np.any(neighbour < -1)
        or np.any(neighbour >= n)
    ):
        raise ValueError("invalid or nonfinite pressure operator/FV topology")
    internal = np.flatnonzero(neighbour >= 0)
    origins, destinations = owner[internal], neighbour[internal]
    if len(internal) == 0 or np.any(origins == destinations):
        raise ValueError("missing or degenerate internal FV faces")
    # Check *both directions* and the complete sparsity graph. This prevents
    # inventing a physical face based on a CSR edge that is not in the mesh.
    directed = np.concatenate((origins * n + destinations, destinations * n + origins))
    if len(np.unique(directed)) != len(directed):
        raise ValueError("multiple physical faces for one matrix adjacency")
    coo = a.tocoo()
    off = coo.row != coo.col
    graph = coo.row[off].astype(np.int64) * n + coo.col[off]
    if len(graph) != len(directed) or not np.array_equal(
        np.sort(graph), np.sort(directed)
    ):
        raise ValueError("CSR off-diagonal graph differs from real FV faces")
    forward = np.asarray(a[origins, destinations]).ravel()
    reverse = np.asarray(a[destinations, origins]).ravel()
    if (
        np.any(forward >= 0.0)
        or np.any(reverse >= 0.0)
        or not np.allclose(forward, reverse, rtol=1e-12, atol=0.0)
    ):
        raise ValueError("pressure CSR is not a symmetric M-matrix on FV faces")
    dtype = np.longdouble
    row_sums = np.empty(n, dtype=dtype)
    row_abs = np.empty(n, dtype=dtype)
    for i in range(n):
        values = a.data[a.indptr[i] : a.indptr[i + 1]].astype(dtype)
        row_sums[i] = np.sum(values, dtype=dtype)
        row_abs[i] = np.sum(np.abs(values), dtype=dtype)
    anchored = row_sums > row_abs * dtype(1e-11)
    if np.any(row_sums < -row_abs * dtype(1e-11)):
        raise ValueError("pressure matrix has negative boundary leakage")
    if (
        len(outlet_faces) != int(np.count_nonzero(anchored))
        or np.any(neighbour[outlet_faces] != -1)
        or len(np.unique(outlet_faces)) != len(outlet_faces)
    ):
        raise ValueError("pressure anchor count != actual Dirichlet outlet faces")
    boundary_owner = owner[outlet_faces]
    if len(np.unique(boundary_owner)) != len(boundary_owner) or not np.array_equal(
        np.sort(boundary_owner), np.flatnonzero(anchored)
    ):
        raise ValueError("pressure anchor cells differ from outlet face owners")
    h = high.astype(dtype)
    l = low.astype(dtype)
    # True physical-face owner-to-neighbour orientation and original matrix
    # coefficient. The face response is one number, shared with opposite signs
    # by owner and neighbour. The nonorthogonal Rhie-Chow flux is NOT here.
    face_internal = -forward.astype(dtype) * (
        (h[origins] - h[destinations]) + (l[origins] - l[destinations])
    )
    face_boundary = row_sums[boundary_owner] * (h[boundary_owner] + l[boundary_owner])
    divergence = np.zeros(n, dtype=dtype)
    np.add.at(divergence, origins, face_internal)
    np.add.at(divergence, destinations, -face_internal)
    internal_global_sum = np.sum(divergence, dtype=dtype)
    np.add.at(divergence, boundary_owner, face_boundary)
    residual, _, eligible = compensated_residual(a, rhs, high, low)
    if not eligible:
        raise ValueError("compensated original CSR is not symmetric positive-diagonal")
    scale = np.sqrt(np.sum(rhs.astype(dtype) ** 2, dtype=dtype))
    if not scale > 0:
        raise ValueError("zero RHS cannot support a relative-residual certificate")
    face_defect = rhs.astype(dtype) - divergence
    original = np.sqrt(np.sum(residual**2, dtype=dtype)) / scale
    face_norm = np.sqrt(np.sum(face_defect**2, dtype=dtype)) / scale
    identity = np.max(np.abs(face_defect - residual))
    # Do not equate linear balance with complete physical pressure correction:
    # one cannot infer the original predictor flux from an A*x=b capture.
    return {
        "schema": "cfdx-bfs-implicit-fv-pressure-face-audit/v1",
        "status": "IMPLICIT_PRESSURE_FACES_ONLY_NOT_FULL_RHIE_CHOW",
        "cells": n,
        "internal_fv_faces": len(internal),
        "outlet_pressure_faces": len(outlet_faces),
        "csr_nonzeros": int(a.nnz),
        "csr_face_topology_bijection": True,
        "relative_original_compensated_residual": float(original),
        "relative_real_mesh_implicit_face_defect": float(face_norm),
        "max_cell_operator_face_identity": float(identity),
        "internal_face_pair_net": float(internal_global_sum),
        "unchanged_original_tolerance": float(target_tolerance),
        "original_csr_linear_gate_pass": bool(original <= target_tolerance),
        "real_mesh_implicit_face_gate_pass": bool(face_norm <= target_tolerance),
        "implicit_face_linear_gate_pass": bool(
            original <= target_tolerance and face_norm <= target_tolerance
        ),
        "actual_rhie_chow_predictor_and_nonorth_flux_replayed": False,
        "cell_velocity_and_pressure_update_replayed": False,
        "mesh_capture_cryptographic_binding_proven": False,
        "physical_fv_face_flux_certified": False,
        "sst_wall_shear_reattachment_qualified": False,
    }


def refine_implicit_faces(
    matrix: csr_matrix,
    rhs: np.ndarray,
    high: np.ndarray,
    low: np.ndarray,
    owner: np.ndarray,
    neighbour: np.ndarray,
    outlet_faces: np.ndarray,
    *,
    target_tolerance: float,
    iterations: int,
) -> dict:
    """Offline LU refinement of REAL FV topology's implicit pressure action.

    This is NOT the CFDX production solver. It deliberately tracks both the
    original rounded CSR residual and the physical-face implicit residual:
    improving one must never be represented as satisfying the other.
    """
    if not 0 <= iterations <= 3:
        raise ValueError("FV face diagnostic refinement budget must be 0..3")
    high, low = high.copy(), low.copy()
    before = verify_face_operator(
        matrix,
        rhs,
        high,
        low,
        owner,
        neighbour,
        outlet_faces,
        target_tolerance=target_tolerance,
    )
    observations = [before]
    if iterations:
        # A double-precision LU supplies only the directional correction;
        # accurate accumulation and the face residual remain high/low.
        factor = splu(matrix.tocsc())
        abs_rows = np.asarray(abs(matrix).sum(axis=1)).ravel().astype(np.longdouble)
        for _ in range(iterations):
            rcsr, rows, _ = compensated_residual(matrix, rhs, high, low)
            # For a non-Dirichlet cell the O(eps) row sum of the rounded
            # diagonal is NOT a boundary flux. Subtract its spurious action
            # from the CSR operator to match the real owner/neighbour faces.
            anchored = rows > abs_rows * np.longdouble(1e-11)
            face_residual = rcsr + np.where(
                ~anchored,
                rows * (high.astype(np.longdouble) + low.astype(np.longdouble)),
                np.longdouble(0.0),
            )
            delta = factor.solve(np.asarray(face_residual, dtype=np.float64))
            if not np.isfinite(delta).all():
                raise ArithmeticError("nonfinite FV pressure refinement")
            two_sum_update(high, low, delta)
            now = verify_face_operator(
                matrix,
                rhs,
                high,
                low,
                owner,
                neighbour,
                outlet_faces,
                target_tolerance=target_tolerance,
            )
            observations.append(now)
            if now["relative_real_mesh_implicit_face_defect"] > (
                observations[-2]["relative_real_mesh_implicit_face_defect"] * 1.01
            ):
                raise ArithmeticError(
                    "FV pressure refinement increased the true face defect"
                )
    final = observations[-1]
    return {
        "schema": "cfdx-bfs-implicit-face-refinement/v1",
        "status": "OFFLINE_IMPLICIT_FV_FACE_DIAGNOSTIC_ONLY",
        "refinement_count": iterations,
        "observations": observations,
        "final_implicit_face_gate_pass": final["implicit_face_linear_gate_pass"],
        "final_real_mesh_implicit_face_gate_pass": final[
            "real_mesh_implicit_face_gate_pass"
        ],
        "original_csr_and_face_operator_conflict": bool(
            final["real_mesh_implicit_face_gate_pass"]
            and not final["original_csr_linear_gate_pass"]
        ),
        "final_original_csr_residual": final["relative_original_compensated_residual"],
        "no_production_state_or_rhie_chow_flux_modified": True,
        "physical_fv_flux_qualified": False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--split", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pressure-patch", default="outlet")
    parser.add_argument(
        "--refinements",
        type=int,
        default=0,
        help="offline face-residual LU corrections (0..3), never production",
    )
    args = parser.parse_args(argv)
    a, rhs, meta = read_capture(args.capture)
    checksum = args.capture.read_text(encoding="ascii").splitlines()[-1].strip()
    verify_capture_checksum(args.capture, a, rhs, checksum)
    high, low = read_split(args.split, a.shape[0], checksum)
    n, owner, neighbour, outlet = read_mesh(args.mesh, args.pressure_patch)
    if n != a.shape[0]:
        raise ValueError("mesh and captured pressure matrix have different sizes")
    report = refine_implicit_faces(
        a,
        rhs,
        high,
        low,
        owner,
        neighbour,
        outlet,
        target_tolerance=meta["target_tolerance"],
        iterations=args.refinements,
    )
    report["provenance"] = {
        "mesh_sha256": sha256(args.mesh),
        "capture_sha256": sha256(args.capture),
        "split_sha256": sha256(args.split),
        "split_capture_checksum": checksum,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["final_implicit_face_gate_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
