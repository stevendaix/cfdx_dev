#!/usr/bin/env python3
"""Validate controlled wall-distance benchmark results.

This is a numerical qualification gate: execution alone is not a PASS.
Geometric/graph methods are references or propagation diagnostics and are
therefore excluded from the PDE convergence gate.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict

PDE_METHODS = {
    "poisson": 1e-6,
    "eikonal": 1e-6,
    "hamilton_jacobi": 1e-6,
    "advection_diffusion": 1e-6,
    "hybrid_poisson_eikonal": 1e-6,
}

def main(path: str) -> int:
    errors: list[str] = []
    rows_by_method: dict[str, list[dict[str, str]]] = defaultdict(list)

    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    required = {
        "N", "method", "residual_inf", "monotonicity_violations",
        "iterations", "converged", "stopping_reason",
    }
    if not rows:
        print("ERROR: controlled benchmark CSV is empty")
        return 1
    missing = required - set(rows[0])
    if missing:
        print("ERROR: controlled benchmark CSV missing columns:", ", ".join(sorted(missing)))
        return 1

    for row in rows:
        method = row["method"]
        if method in PDE_METHODS:
            rows_by_method[method].append(row)

    for method, tolerance in PDE_METHODS.items():
        method_rows = rows_by_method.get(method, [])
        if not method_rows:
            errors.append(f"[{method}] no benchmark rows found")
            continue

        for row in method_rows:
            n = int(row["N"])
            residual = float(row["residual_inf"])
            violations = float(row["monotonicity_violations"])
            iterations = int(row["iterations"])
            converged = row["converged"].strip().lower() == "true"
            reason = row["stopping_reason"]

            if not converged:
                errors.append(
                    f"[{method}] N={n}: not converged "
                    f"(iterations={iterations}, reason={reason}, residual={residual:.6e})"
                )
            if not residual == residual or residual > tolerance:
                errors.append(
                    f"[{method}] N={n}: residual {residual:.6e} > {tolerance:.6e}"
                )

            # The metric counts violated directed neighbour relations. The
            # controlled benchmark must expose essentially monotone distance.
            # Use a strict relative gate, but avoid treating tiny integer
            # roundoff counts as a failure on coarse grids.
            samples = max(1, n ** 3)
            if violations > 0.01 * samples:
                errors.append(
                    f"[{method}] N={n}: monotonicity violations={violations:.0f} "
                    f"exceed 1% of grid-cell count"
                )

    if errors:
        print("CONTROLLED WALL-DISTANCE NUMERICAL GATE: FAIL")
        for error in errors:
            print(" -", error)
        return 1

    print("CONTROLLED WALL-DISTANCE NUMERICAL GATE: PASS")
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
