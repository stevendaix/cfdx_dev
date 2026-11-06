"""Fail-closed early-iteration NACA physical-history equivalence gate.

Compare two production outer histories at the SAME configured pressure solver,
mesh, iteration count and initial state. This is not a final force validation.
"""

import argparse
import csv
import math
from pathlib import Path
import sys

PHYSICS = (
    "continuity_normalized", "momentum_rhs_relative",
    "mass_normalized_imbalance", "mass_local_linf",
    "flux_velocity_mismatch_normalized",
)
COUNTS = (
    "iteration", "pressure_linear_iterations", "momentum_linear_iterations",
    "nonfinite_mass_faces", "nonfinite_velocities",
)


def compare_histories(base: str, candidate: str, max_abs: float = 1e-12,
                      expected_iterations: int = 3) -> dict:
    """Compare actual output CSV contents, never assume a matching solver config."""
    from io import StringIO
    before = list(csv.DictReader(StringIO(base)))
    after = list(csv.DictReader(StringIO(candidate)))
    if len(before) != expected_iterations or len(after) != expected_iterations:
        raise ValueError("different or incomplete outer histories")
    maximum = {metric: 0.0 for metric in PHYSICS}
    for count, (row_before, row_after) in enumerate(zip(before, after), start=1):
        for name in COUNTS:
            if int(row_before[name]) != int(row_after[name]):
                raise ValueError(f"{name} differs on outer iteration {count}")
        if int(row_before["iteration"]) != count:
            raise ValueError("unexpected or reordered outer iterations")
        for name in PHYSICS:
            a, b = float(row_before[name]), float(row_after[name])
            if not math.isfinite(a) or not math.isfinite(b):
                raise ValueError(f"nonfinite metric {name}")
            difference = abs(a - b)
            maximum[name] = max(maximum[name], difference)
            if difference > max_abs:
                raise ValueError(f"{name} differs by {difference:g} on iteration {count}")
    return maximum


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--expected-iterations", type=int, default=3)
    parser.add_argument("--max-absolute", type=float, default=1e-12)
    args = parser.parse_args()
    if args.expected_iterations < 1 or not math.isfinite(args.max_absolute) or args.max_absolute < 0:
        parser.error("invalid iteration count or absolute tolerance")
    maximum = compare_histories(args.baseline.read_text(encoding="utf-8"),
                                args.candidate.read_text(encoding="utf-8"),
                                max_abs=args.max_absolute,
                                expected_iterations=args.expected_iterations)
    for name, value in maximum.items():
        print(f"{name}: max_abs_delta={value:.17g}")
    print("PASS: matching NACA outer histories on configured early iterations")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, OverflowError) as exc:
        print("NACA_PHYSICAL_HISTORY_FAIL:", exc, file=sys.stderr)
        raise SystemExit(1) from exc
