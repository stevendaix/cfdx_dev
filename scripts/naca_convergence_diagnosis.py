#!/usr/bin/env python3
"""Strict per-component NACA history diagnostics; NEVER a convergence certificate.

A valid history is not evidence that the corresponding solver returned success.
The official converged flag, complete force certificate, bounded high-order
scheme and 3-grid results are separate mandatory qualification gates.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

REQUIRED = (
    "iteration",
    "momentum_rhs_relative",
    "momentum_reference_relative",
    "momentum_x_rhs_relative",
    "momentum_y_rhs_relative",
    "momentum_z_rhs_relative",
    "momentum_component_max",
    "continuity_normalized",
    "face_momentum_x_relative",
    "face_momentum_y_relative",
    "face_momentum_z_relative",
    "velocity_change_inf",
    "pressure_change_inf",
    "pressure_residual",
    "mass_normalized_imbalance",
    "high_order_fraction",
    "convection_scheme",
    "algorithm",
)
NONNEGATIVE = tuple(
    k for k in REQUIRED if k not in ("iteration", "convection_scheme", "algorithm")
)


def read_history(path: Path) -> list[dict[str, float | int | str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        names = tuple(reader.fieldnames or ())
        if len(names) != len(set(names)) or not set(REQUIRED).issubset(names):
            raise ValueError("NACA diagnostic requires all per-component columns")
        rows: list[dict[str, float | int | str]] = []
        for expected, row in enumerate(reader, start=1):
            if (
                row is None
                or None in row
                or any(value is None for value in row.values())
            ):
                raise ValueError("NACA history has malformed CSV row")
            raw_iter = row["iteration"]
            if not raw_iter or not raw_iter.isdecimal() or int(raw_iter) != expected:
                raise ValueError("NACA iteration history is incomplete or out of order")
            item: dict[str, float | int | str] = {"iteration": expected}
            for field in NONNEGATIVE:
                value = float(row[field])
                if not math.isfinite(value) or value < 0.0:
                    raise ValueError(
                        f"NACA history has invalid nonnegative metric: {field}"
                    )
                item[field] = value
            if float(item["high_order_fraction"]) > 1.0:
                raise ValueError("NACA history has invalid high-order fraction")
            for field, options in (
                ("algorithm", ("SIMPLE", "SIMPLEC")),
                (
                    "convection_scheme",
                    ("UPWIND", "TVD", "TVD_VANLEER", "TVD_VANALBADA"),
                ),
            ):
                if row[field] not in options:
                    raise ValueError(f"NACA history has invalid {field}")
                item[field] = row[field]
            components = [
                float(item[f"momentum_{axis}_rhs_relative"]) for axis in "xyz"
            ]
            advertised = float(item["momentum_component_max"])
            if not math.isclose(
                max(components), advertised, abs_tol=1e-10, rel_tol=1e-7
            ):
                raise ValueError(
                    "NACA history momentum component maximum is inconsistent"
                )
            rows.append(item)
    if not rows:
        raise ValueError("NACA convergence history is empty")
    return rows


def diagnose(rows: list[dict[str, float | int | str]], tolerance: float = 1e-8) -> dict:
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("invalid nonlinear tolerance")
    last = rows[-1]
    scalar_checks = (
        "momentum_rhs_relative",
        "momentum_x_rhs_relative",
        "momentum_y_rhs_relative",
        "momentum_z_rhs_relative",
        "face_momentum_x_relative",
        "face_momentum_y_relative",
        "face_momentum_z_relative",
        "continuity_normalized",
        "velocity_change_inf",
        "pressure_change_inf",
        "pressure_residual",
    )
    blockers = {
        key: float(last[key]) for key in scalar_checks if float(last[key]) > tolerance
    }
    # A new restart may legitimately change the reference residual. Its
    # initial-relative ratio is a separate solver policy, NOT an oracle here.
    ratio = float(last["momentum_reference_relative"])
    span = rows[-min(8, len(rows)) :]
    trend = None
    if len(span) >= 2:
        before = float(span[0]["momentum_component_max"])
        after = float(span[-1]["momentum_component_max"])
        if before > 0.0:
            trend = (after / before) ** (1 / (len(span) - 1))
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_A_SOLVER_CERTIFICATE",
        "iterations": int(last["iteration"]),
        "tolerance_unchanged": tolerance,
        "algorithm": last["algorithm"],
        "scheme": last["convection_scheme"],
        "full_order_fraction": float(last["high_order_fraction"]),
        "full_order_scheme": last["convection_scheme"]
        in ("TVD", "TVD_VANLEER", "TVD_VANALBADA")
        and float(last["high_order_fraction"]) == 1.0,
        "momentum_reference_ratio_diagnostic": ratio,
        "momentum_component_max_geometric_ratio": trend,
        "unmet_nonlinear_checks": blockers,
        "physical_tolerance_checks_pass": len(blockers) == 0,
        "solver_convergence_proven": False,
        "three_grid_gci_proven": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = diagnose(read_history(args.history), args.tolerance)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        if args.output.resolve() == args.history.resolve():
            parser.error("diagnostic output must not overwrite the history")
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
