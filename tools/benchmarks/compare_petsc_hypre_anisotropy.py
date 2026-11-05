"""Independent reference for CFDX's native anisotropic pressure benchmarks.

The production CFDX path has NO PETSc/HYPRE runtime dependency. Run this
optional probe in a separate Python environment with petsc4py, PETSc and
HYPRE enabled to compare native Boomer-style AMG against actual PCHYPRE
BoomerAMG and native smoothed aggregation against PETSc PCGAMG.

Without petsc4py, --backend petsc explicitly exits with 77 (SKIP); it never
reports a fictitious validation. --backend scipy uses an independent direct
sparse factorization and is useful even without PETSc/HYPRE installed.

Usage:
 python3 tools/benchmarks/compare_petsc_hypre_anisotropy.py --backend scipy --nx 64 --ny 64
 python3 tools/benchmarks/compare_petsc_hypre_anisotropy.py --backend petsc --nx 64 --ny 64
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time


def matrix_entries(nx: int, ny: int, ratio: float, axis: str):
    ax, ay = (ratio, 1.0) if axis == "x" else (1.0, ratio)
    for j in range(ny):
        for i in range(nx):
            k = j * nx + i
            diagonal = 1.0  # reaction term exactly matches the C++ regression
            neighbors: list[tuple[int, float]] = []
            for valid, neighbor, weight in (
                (i > 0, k - 1, ax),
                (i + 1 < nx, k + 1, ax),
                (j > 0, k - nx, ay),
                (j + 1 < ny, k + nx, ay),
            ):
                if valid:
                    neighbors.append((neighbor, -weight))
                    diagonal += weight
            yield k, neighbors + [(k, diagonal)]


def exact_solution(nx: int, ny: int):
    return [
        0.25 + math.sin(math.pi * (i + 0.5) / nx)
        * math.sin(2.0 * math.pi * (j + 0.5) / ny)
        for j in range(ny)
        for i in range(nx)
    ]


def reference_scipy(nx: int, ny: int, ratio: float, axis: str):
    try:
        import numpy as np
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import spsolve
    except ImportError as error:
        raise RuntimeError("SciPy/numpy missing; independent reference unavailable") from error
    indices = []
    columns = []
    values = []
    for i, entries in matrix_entries(nx, ny, ratio, axis):
        for j, value in entries:
            indices.append(i)
            columns.append(j)
            values.append(value)
    A = csr_matrix((values, (indices, columns)), shape=(nx * ny, nx * ny))
    exact = np.asarray(exact_solution(nx, ny))
    b = A @ exact
    setup_start = time.perf_counter()
    # SciPy's sparse direct solve performs its factorization in solve; do not
    # incorrectly attribute that work to an AMG setup stage.
    setup_seconds = time.perf_counter() - setup_start
    start = time.perf_counter()
    x = spsolve(A, b)
    solve_seconds = time.perf_counter() - start
    relative_residual = float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))
    relative_error = float(np.linalg.norm(x - exact) / np.linalg.norm(exact))
    return "scipy_spsolve", 0, relative_residual, relative_error, setup_seconds, solve_seconds


def reference_petsc(nx: int, ny: int, ratio: float, axis: str, ksp_type: str, pc_type: str):
    from petsc4py import PETSc
    n = nx * ny
    A = PETSc.Mat().createAIJ(size=(n, n), nnz=5, comm=PETSc.COMM_SELF)
    for i, entries in matrix_entries(nx, ny, ratio, axis):
        A.setValues(i, [col for col, _ in entries], [value for _, value in entries])
    A.assemble()
    known = A.createVecRight()
    known.setValues(list(range(n)), exact_solution(nx, ny))
    known.assemble()
    rhs = A.createVecLeft()
    A.mult(known, rhs)
    x = A.createVecRight()
    x.set(0)
    ksp = PETSc.KSP().create(PETSc.COMM_SELF)
    ksp.setOperators(A)
    ksp.setType(ksp_type)
    if ksp_type in ("gmres", "fgmres"):
        ksp.setGMRESRestart(35)  # Same restart as the native CFDX benchmark.
    pc = ksp.getPC()
    if pc_type == "hypre_boomeramg":
        pc.setType(PETSc.PC.Type.HYPRE)
        pc.setHYPREType("boomeramg")
    elif pc_type == "gamg":
        pc.setType(PETSc.PC.Type.GAMG)
    elif pc_type == "ilu":
        pc.setType(PETSc.PC.Type.ILU)
    else:
        raise ValueError(f"unsupported PC: {pc_type}")
    ksp.setTolerances(rtol=1e-9, atol=0.0, max_it=350)
    ksp.setFromOptions()
    try:
        begin = time.perf_counter()
        ksp.setUp()
        setup_seconds = time.perf_counter() - begin
        begin = time.perf_counter()
        ksp.solve(rhs, x)
        solve_seconds = time.perf_counter() - begin
        residual = A.createVecLeft()
        A.mult(x, residual)
        residual.axpy(-1.0, rhs)
        diff = x.duplicate()
        x.copy(diff)
        diff.axpy(-1.0, known)
        relative_residual = residual.norm() / max(rhs.norm(), 1e-300)
        relative_error = diff.norm() / max(known.norm(), 1e-300)
        return (f"petsc_{ksp_type}_{pc_type}", int(ksp.getIterationNumber()),
                relative_residual, relative_error, setup_seconds, solve_seconds,
                int(ksp.getConvergedReason()))
    finally:
        ksp.destroy()
        x.destroy()
        rhs.destroy()
        known.destroy()
        A.destroy()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--backend", choices=("scipy", "petsc"), default="scipy")
    p.add_argument("--nx", type=int, default=64)
    p.add_argument("--ny", type=int, default=64)
    p.add_argument("--ratios", nargs="*", type=float, default=[1, 1e3, 1e6])
    p.add_argument("--axes", nargs="*", choices=("x", "y"), default=["x", "y"])
    args = p.parse_args()
    if args.nx < 4 or args.ny < 4 or args.nx * args.ny > 262144:
        p.error("mesh size outside qualification range")
    if any(not math.isfinite(r) or r < 1 for r in args.ratios):
        p.error("all ratios must be finite and >= 1")
    if args.backend == "petsc":
        try:
            from petsc4py import PETSc  # noqa: F401
        except (ImportError, RuntimeError) as error:
            print(f"SKIP(77): PETSc/HYPRE not installed: {error}", file=sys.stderr)
            return 77

    writer = csv.writer(sys.stdout)
    writer.writerow(["source", "axis", "ratio", "n", "iterations", "true_residual",
                     "relative_solution_error", "setup_s", "solve_s", "reason"])
    attempted = qualified = 0
    for axis in args.axes:
        for ratio in args.ratios:
            configurations = (
                [("scipy", "direct")]
                if args.backend == "scipy" else
                [("fgmres", "hypre_boomeramg"), ("fgmres", "gamg"),
                 ("gmres", "ilu"), ("fgmres", "ilu"), ("cg", "gamg")]
            )
            for krylov, preconditioner in configurations:
                try:
                    attempted += 1
                    if args.backend == "scipy":
                        name, iterations, residual, err, setup, solve = reference_scipy(
                            args.nx, args.ny, ratio, axis)
                        reason = "reference"
                    else:
                        (name, iterations, residual, err, setup, solve,
                         convergence_reason) = reference_petsc(
                            args.nx, args.ny, ratio, axis, krylov, preconditioner)
                        reason = str(convergence_reason)
                    good = (math.isfinite(residual) and math.isfinite(err)
                            and residual <= 1e-9 and err <= 2e-3
                            and (args.backend == "scipy" or int(reason) > 0))
                    if good:
                        qualified += 1
                    writer.writerow([name, axis, ratio, args.nx * args.ny, iterations,
                                     f"{residual:.12g}", f"{err:.12g}",
                                     f"{setup:.6g}", f"{solve:.6g}", reason])
                except Exception as error:
                    # A missing HYPRE PC in a PETSc build is a SKIP for this
                    # configuration, not a positive qualification certificate.
                    print(f"NOT VALIDATED {args.backend}/{krylov}/{preconditioner} "
                          f"axis={axis} ratio={ratio}: {error}", file=sys.stderr)
    print(f"qualification: {qualified}/{attempted} reference combinations passed", file=sys.stderr)
    if qualified == 0:
        return 77
    return 0 if qualified == attempted else 1


if __name__ == "__main__":
    raise SystemExit(main())