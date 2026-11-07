#!/usr/bin/env python3
"""Localize NACA momentum defects without treating linear convergence as CFD PASS.

Consumes the accepted-iteration CSV exported by test_naca0012_qualification.
The report is always a DIAGNOSTIC; only the CFD solver's independent full
nonlinear convergence and mesh-qualification gates can certify the airfoil.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

SCALARS = (
    "momentum_rhs_relative", "continuity_normalized", "mass_local_linf",
    "momentum_rhs_x", "momentum_rhs_y", "momentum_rhs_z",
    "momentum_balance_x", "momentum_balance_y", "momentum_balance_z",
    "alpha_u", "alpha_p", "velocity_change_inf", "pressure_change_inf",
)
COUNTS = ("momentum_worst_cell_x", "momentum_worst_cell_y",
          "momentum_worst_cell_z", "nonfinite_mass_faces", "nonfinite_velocities",
          "pressure_linear_iterations", "momentum_linear_iterations")


def read_history(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf8") as f:
        reader = csv.DictReader(f)
        required = set(SCALARS + COUNTS + ("iteration",))
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"incomplete accepted-iteration CSV fields: {sorted(required-set(reader.fieldnames or []))}")
        values = []
        last_it = 0
        for raw in reader:
            row = {k: float(raw[k]) for k in SCALARS}
            row.update({k: int(raw[k]) for k in COUNTS})
            row["iteration"] = int(raw["iteration"])
            if row["iteration"] <= last_it:
                raise ValueError("iteration numbers are not strictly increasing")
            last_it = row["iteration"]
            if any(not math.isfinite(row[k]) for k in SCALARS):
                raise ValueError(f"nonfinite momentum/continuity at iteration {last_it}")
            if any(row[k] < 0 for k in COUNTS):
                raise ValueError("negative iteration count, cell or nonfinite counter")
            if row["nonfinite_mass_faces"] or row["nonfinite_velocities"]:
                raise ValueError(f"nonfinite solver state at accepted iteration {last_it}")
            if (row["alpha_u"] <= 0 or row["alpha_u"] > 1 or
                    row["alpha_p"] <= 0 or row["alpha_p"] > 1):
                raise ValueError("invalid effective relaxation factor")
            values.append(row)
    if not values:
        raise ValueError("empty accepted-iteration history")
    return values


def diagnose(rows: list[dict], physical_gate: float = 1e-8) -> dict:
    if not math.isfinite(physical_gate) or physical_gate <= 0:
        raise ValueError("invalid physical threshold")
    # On a planar mesh, z-component is diagnostically reported but is NOT a
    # physical 3-D transport equation to be used for artificial stabilization.
    worst = max(rows, key=lambda r: max(r["momentum_balance_x"], r["momentum_balance_y"]))
    axis = "x" if worst["momentum_balance_x"] >= worst["momentum_balance_y"] else "y"
    final = rows[-1]
    pressure_leads = next((r["iteration"] for r in rows
        if r["continuity_normalized"] <= physical_gate and
           r["momentum_rhs_relative"] > physical_gate), None)
    window = rows[-min(10, len(rows)):]
    baseline = window[0]["momentum_rhs_relative"]
    last = window[-1]["momentum_rhs_relative"]
    # Do not label a plateau from one or two iterations alone.
    stalled = len(window) >= 6 and baseline > physical_gate and last > 0 and last/baseline > 0.95
    rhs_scale_gap = max(final["momentum_rhs_x"], final["momentum_rhs_y"]) / max(
        final["momentum_rhs_relative"], 1e-300)
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_PHYSICALLY_QUALIFIED",
        "accepted_iterations": len(rows),
        "last_iteration": final["iteration"],
        "final_momentum_rhs_relative": final["momentum_rhs_relative"],
        "final_continuity_normalized": final["continuity_normalized"],
        "final_local_mass_linf": final["mass_local_linf"],
        "first_pressure_ahead_of_momentum_iteration": pressure_leads,
        "momentum_plateau_last_window": stalled,
        "strong_component_rhs_scale_disparity": rhs_scale_gap > 100.0,
        "component_to_global_rhs_ratio": rhs_scale_gap,
        "worst_independent_momentum_balance": {
            "iteration": worst["iteration"], "component": axis,
            "cell": worst[f"momentum_worst_cell_{axis}"],
            "normalized_defect": worst[f"momentum_balance_{axis}"],
        },
        "suggested_next_diagnostic": (
            "inspect face traction, boundary gradient and Rhie-Chow flux at "
            "the indicated cell; compare local residual before and after correction"
        ),
        "note": "A pressure Krylov PASS cannot promote a non-converged momentum solution",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("history", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    try:
        report = diagnose(read_history(args.history))
        temp = args.output.with_name(args.output.name + ".tmp")
        temp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")
        os.replace(temp, args.output)
    except (OSError, ValueError, KeyError, OverflowError) as exc:
        args.output.unlink(missing_ok=True)
        print(f"NACA_NONLINEAR_DIAGNOSTIC: FAIL {exc}")
        return 2
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
