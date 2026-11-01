"""Estimate the Jacobi-scaled pressure-matrix spectrum from a Matrix-Market dump.

This is a forensic diagnostic, not a qualification gate. It checks matrix
symmetry, graph connectivity and positive diagonal before estimating the
extreme eigenvalues of D^{-1/2} A D^{-1/2} with ARPACK.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any


def parse_rhs_csv(path: Path) -> list[dict[str, float | int]]:
    """Read the row-to-physical-cell mapping emitted by CFDX."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"row", "cell", "rhs"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: expected CSV columns {sorted(required)}")
        records = []
        for record in reader:
            records.append(
                {
                    "row": int(record["row"]),
                    "cell": int(record["cell"]),
                    "rhs": float(record["rhs"]),
                }
            )
    return records


def condition_estimate(lambda_min: float, lambda_max: float) -> float | None:
    """Return a condition estimate only for a finite positive spectrum."""
    if (
        not math.isfinite(lambda_min)
        or not math.isfinite(lambda_max)
        or lambda_min <= 0.0
        or lambda_max <= 0.0
    ):
        return None
    return lambda_max / lambda_min


def _matrix_inf_norm(matrix: Any) -> float:
    import numpy as np

    return float(np.asarray(abs(matrix).sum(axis=1)).ravel().max(initial=0.0))


def analyze_matrix(path: Path, rhs_path: Path | None = None) -> dict[str, Any]:
    """Analyze a sparse pressure matrix and return diagnostic quantities."""
    try:
        import numpy as np
        from scipy.io import mmread
        from scipy.sparse import csr_matrix, diags, issparse
        from scipy.sparse.csgraph import connected_components
        from scipy.sparse.linalg import ArpackNoConvergence, eigsh
    except ImportError as exc:
        raise RuntimeError(
            "spectrum analysis requires NumPy and SciPy (scipy.sparse/eigsh)"
        ) from exc

    matrix = mmread(path)
    matrix = matrix.tocsr() if issparse(matrix) else csr_matrix(matrix)
    if matrix.shape[0] != matrix.shape[1] or matrix.shape[0] == 0:
        raise ValueError(
            f"{path}: expected a non-empty square matrix, got {matrix.shape}"
        )
    if not np.isfinite(matrix.data).all():
        raise ValueError(f"{path}: matrix contains non-finite coefficients")

    n = matrix.shape[0]
    diagonal = np.asarray(matrix.diagonal(), dtype=float)
    if not np.isfinite(diagonal).all() or np.any(diagonal <= 0.0):
        raise ValueError("Jacobi scaling requires a finite, strictly positive diagonal")

    norm_inf = _matrix_inf_norm(matrix)
    symmetry_defect = _matrix_inf_norm(matrix - matrix.T)
    symmetry_relative = symmetry_defect / max(norm_inf, np.finfo(float).tiny)

    graph = matrix.copy()
    graph.setdiag(0.0)
    graph.eliminate_zeros()
    graph.data[:] = 1.0
    graph = graph.maximum(graph.T)
    component_count, labels = connected_components(
        graph, directed=False, return_labels=True
    )

    inv_sqrt_diag = 1.0 / np.sqrt(diagonal)
    scale = diags(inv_sqrt_diag)
    scaled = (scale @ matrix @ scale).tocsr()
    scaled_symmetry = _matrix_inf_norm(scaled - scaled.T) / max(
        _matrix_inf_norm(scaled), np.finfo(float).tiny
    )

    result: dict[str, Any] = {
        "rows": int(n),
        "nnz": int(matrix.nnz),
        "matrix_inf_norm": norm_inf,
        "symmetry_defect_inf": symmetry_defect,
        "symmetry_relative": symmetry_relative,
        "scaled_symmetry_relative": scaled_symmetry,
        "diagonal_min": float(diagonal.min()),
        "diagonal_max": float(diagonal.max()),
        "connected_components": int(component_count),
        "largest_component_cells": int(np.bincount(labels).max()),
        "lambda_min": None,
        "lambda_max": None,
        "lambda_min_residual": None,
        "lambda_max_residual": None,
        "condition_estimate": None,
    }

    if rhs_path is not None:
        rhs_records = parse_rhs_csv(rhs_path)
        if len(rhs_records) != n:
            raise ValueError(
                f"RHS row count {len(rhs_records)} does not match matrix size {n}"
            )
        rhs_abs = [abs(float(item["rhs"])) for item in rhs_records]
        max_index = max(range(n), key=rhs_abs.__getitem__)
        result["rhs_linf"] = rhs_abs[max_index]
        result["rhs_row"] = int(rhs_records[max_index]["row"])
        result["rhs_cell"] = int(rhs_records[max_index]["cell"])

    # A nonsymmetric or disconnected operator is reported, not silently
    # coerced into a symmetric eigenproblem.
    if scaled_symmetry > 1.0e-10 or component_count != 1:
        result["spectrum_status"] = "NOT_ESTIMATED_MATRIX_STRUCTURE"
        return result

    if n <= 2:
        dense_eigenvalues = np.linalg.eigvalsh(scaled.toarray())
        eigen_min = float(dense_eigenvalues[0])
        eigen_max = float(dense_eigenvalues[-1])
        residual_min = residual_max = 0.0
    else:
        ncv = min(80, n - 1)
        maxiter = max(2000, 20 * n)

        def extreme(which: str) -> tuple[float, float]:
            try:
                values, vectors = eigsh(
                    scaled,
                    k=1,
                    which=which,
                    tol=1.0e-9,
                    maxiter=maxiter,
                    ncv=ncv,
                )
            except ArpackNoConvergence as exc:
                raise RuntimeError(
                    f"ARPACK did not converge for which={which}; "
                    f"partial_eigenvalues={exc.eigenvalues}"
                ) from exc
            eigenvalue = float(values[0])
            vector = vectors[:, 0]
            residual = float(
                np.linalg.norm(scaled @ vector - eigenvalue * vector)
                / max(_matrix_inf_norm(scaled), np.finfo(float).tiny)
            )
            return eigenvalue, residual

        eigen_min, residual_min = extreme("SA")
        eigen_max, residual_max = extreme("LA")

    result.update(
        {
            "lambda_min": eigen_min,
            "lambda_max": eigen_max,
            "lambda_min_residual": residual_min,
            "lambda_max_residual": residual_max,
            "condition_estimate": condition_estimate(eigen_min, eigen_max),
            "spectrum_status": (
                "ESTIMATED"
                if eigen_min > 0.0 and max(residual_min, residual_max) <= 1.0e-6
                else "SPECTRUM_NOT_POSITIVE_OR_RESIDUAL_TOO_LARGE"
            ),
        }
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path, help="CFDX Matrix-Market .mtx file")
    parser.add_argument("--rhs", type=Path, help="matching .rhs.csv file")
    args = parser.parse_args()
    try:
        result = analyze_matrix(args.matrix, args.rhs)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"error: {exc}\n")
    for key, value in result.items():
        print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
