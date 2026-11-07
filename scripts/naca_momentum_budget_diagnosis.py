#!/usr/bin/env python3
"""Fail-closed cellwise NACA momentum budget audit (never a CFD certificate)."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

COLUMNS = (
    "iteration", "component", "rank", "cell", "x", "y", "z",
    "convection", "viscous", "transpose_viscous", "pressure", "body_force",
    "bounded_divergence", "mass_divergence", "collocated_difference",
    "physical_residual", "bounded_operator_residual", "worst_residual",
)
TERMS = (
    "convection", "viscous", "transpose_viscous", "pressure", "body_force",
    "bounded_divergence", "mass_divergence", "collocated_difference",
    "physical_residual", "bounded_operator_residual", "worst_residual",
)


def _close(left: float, right: float) -> bool:
    # Machine-level arithmetic consistency, NOT a physics tolerance gate.
    return abs(left - right) <= 2e-10 * max(1.0, abs(left), abs(right))


def read_budget(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("missing or duplicated cell-budget CSV headers")
        if not set(COLUMNS).issubset(reader.fieldnames):
            raise ValueError("missing required cell-budget columns")
        rows: list[dict] = []
        for line, data in enumerate(reader, 2):
            if None in data or any(v is None or not v for v in data.values()):
                raise ValueError(f"malformed cell-budget CSV row {line}")
            row: dict = {}
            for key in ("iteration", "component", "rank", "cell"):
                if not data[key].isdecimal():
                    raise ValueError(f"invalid integer {key} at line {line}")
                row[key] = int(data[key])
            if row["iteration"] < 1 or row["component"] not in (0, 1, 2):
                raise ValueError(f"invalid iteration or component at line {line}")
            for key in ("x", "y", "z", *TERMS):
                val = float(data[key])
                if not math.isfinite(val):
                    raise ValueError(f"nonfinite {key} at line {line}")
                row[key] = val
            if row["worst_residual"] < 0:
                raise ValueError(f"negative worst residual at line {line}")
            conv = row["convection"] + row["viscous"] + row["transpose_viscous"]
            phys = conv + row["pressure"] + row["body_force"]
            if not _close(phys, row["physical_residual"]):
                raise ValueError(f"cell transport/source sum inconsistent at line {line}")
            if not _close(phys + row["bounded_divergence"],
                          row["bounded_operator_residual"]):
                raise ValueError(f"bounded operator sum inconsistent at line {line}")
            rows.append(row)
    if not rows:
        raise ValueError("empty cell-budget file")
    groups: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        groups[row["iteration"], row["component"]].append(row)
    for key, group in groups.items():
        if set(r["rank"] for r in group) != set(range(len(group))):
            raise ValueError(f"noncontiguous ranks for iteration/component {key}")
        if len({r["cell"] for r in group}) != len(group):
            raise ValueError(f"duplicate cells for iteration/component {key}")
        ordered = sorted(group, key=lambda x: x["rank"])
        if any(abs(a["physical_residual"]) + 1e-12 < abs(b["physical_residual"])
               for a, b in zip(ordered, ordered[1:])):
            raise ValueError(f"ranks not ordered by local residual for {key}")
        if not _close(abs(ordered[0]["physical_residual"]), ordered[0]["worst_residual"]):
            raise ValueError(f"top residual disagrees with reported maximum for {key}")
        if any(not _close(r["worst_residual"], ordered[0]["worst_residual"])
               for r in ordered):
            raise ValueError(f"inconsistent worst residual within {key}")
    return rows


def diagnose(rows: list[dict]) -> dict:
    groups: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        groups[row["iteration"], row["component"]].append(row)
    latest = max(key[0] for key in groups)
    summary = []
    for iteration, component in sorted(groups):
        top = min(groups[iteration, component], key=lambda row: row["rank"])
        summary.append({
            "iteration": iteration, "component": "xyz"[component],
            "worst_cell": top["cell"], "x": top["x"], "y": top["y"],
            "convection": top["convection"], "viscous": top["viscous"],
            "pressure": top["pressure"], "collocated_flux_difference":
            top["collocated_difference"],
            "physical_residual": top["physical_residual"],
            "bounded_operator_residual": top["bounded_operator_residual"],
        })
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_A_SOLVER_CERTIFICATE",
        "latest_sampled_iteration": latest,
        "sampled_cell_rows": len(rows),
        "budget_series": summary,
        "solver_convergence_proven": False,
        "three_grid_gci_proven": False,
        "sst_thermal_radiation_qualification_proven": False,
        "flux_note": "collocated difference = authoritative corrected flux minus directly interpolated velocity flux; not pure Rhie-Chow",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("budget", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = diagnose(read_budget(args.budget))
    serialized = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        if args.output.resolve() == args.budget.resolve():
            parser.error("cannot overwrite source CSV")
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
