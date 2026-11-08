"""Compare cold and warm NACA runs without ever certifying convergence.

Both histories must use the same mesh file contents, algorithm, and convection
operator. The original solver's full physical/force certificates, an actual
independent second converged solution, and the GCI are still required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from naca_convergence_diagnosis import read_history

METRICS = (
    "momentum_rhs_relative",
    "momentum_component_max",
    "continuity_normalized",
    "face_momentum_y_relative",
    "mass_normalized_imbalance",
    "velocity_change_inf",
    "pressure_change_inf",
    "pressure_residual",
)


def mesh_sha256(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def compare(
    cold_history: Path, warm_history: Path, cold_mesh: Path, warm_mesh: Path
) -> dict:
    cold_hash = mesh_sha256(cold_mesh)
    warm_hash = mesh_sha256(warm_mesh)
    if cold_hash != warm_hash:
        raise ValueError(
            "cold and warm histories must correspond to identical mesh bytes"
        )
    cold, warm = read_history(cold_history), read_history(warm_history)
    n = min(len(cold), len(warm))
    for iteration in range(n):
        a, b = cold[iteration], warm[iteration]
        for key in ("convection_scheme", "algorithm"):
            if a[key] != b[key]:
                raise ValueError(f"different {key} at iteration {iteration + 1}")
        if not math.isclose(
            float(a["high_order_fraction"]),
            float(b["high_order_fraction"]),
            abs_tol=1e-14,
        ):
            raise ValueError("different convection continuation schedules")
    a, b = cold[n - 1], warm[n - 1]
    results = {}
    for name in METRICS:
        x, y = float(a[name]), float(b[name])
        results[name] = {
            "cold": x,
            "warm": y,
            "cold_over_warm": x / y if y > 0.0 else None,
        }
    return {
        "status": "DIAGNOSTIC_ONLY_NEVER_A_CONVERGENCE_CERTIFICATE",
        "seed_independence_proven": False,
        "mesh_sha256": cold_hash,
        "common_iterations": n,
        "cold_history_length": len(cold),
        "warm_history_length": len(warm),
        "convection_scheme": a["convection_scheme"],
        "algorithm": a["algorithm"],
        "metrics": results,
        "limitations": [
            (
                Mesh equivalence is checked by hashing the two declared inputs; "
                "the CSV format does not authenticate which mesh created each history."
            ),
            "Different RHS normalization scales can affect per-component comparisons.",
            (
                Only independently converged cold and warm solutions, checked for "
                "field and force agreement, can establish initialization independence."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cold-history", type=Path, required=True)
    parser.add_argument("--warm-history", type=Path, required=True)
    parser.add_argument("--cold-mesh", type=Path, required=True)
    parser.add_argument("--warm-mesh", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        result = compare(
            args.cold_history, args.warm_history, args.cold_mesh, args.warm_mesh
        )
        outputs = (
            args.cold_history.resolve(),
            args.warm_history.resolve(),
            args.cold_mesh.resolve(),
            args.warm_mesh.resolve(),
        )
        if args.report and args.report.resolve() in outputs:
            raise ValueError("report must not overwrite an input")
        data = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.report:
            args.report.write_text(data, encoding="utf-8")
        else:
            print(data, end="")
    except (OSError, ValueError, TypeError, OverflowError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
