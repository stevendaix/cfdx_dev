"""Fail-closed V33 SST row forensics. A diagnostic, NEVER a convergence certificate.

Sample: CFDX_SST_ROW_TRACE_V33=341,329 CFDX_SST_DIAGNOSTICS=1 <solver>
  python scripts/analyze_sst_rows_v33.py solver.log --cells 341,329 --output report.json

A valid report proves row trace completeness only. It cannot qualify OpenFOAM,
BFS, nonlinear convergence or one of the production Schur variants.
"""
import argparse
import json
import math
import re
from pathlib import Path

PREFIX = "SST_ROW_TRACE_V33 "
REQUIRED = {"iteration", "cell", "field", "full_residual", "closure_error",
            "source_identity_error", "qualified", "production_limited",
            "omega_viscosity_branch", "branch_changed", "local_relative",
            "winner_accepted", "frozen_source_derivative_valid",
            "dRk_dlogk_source", "dRk_dlogomega_source",
            "dRk_dlogk_fv_source", "dRk_dlogomega_fv_source",
            "source_units", "cell_volume"}


def parse_log(content: str, cells: list[int], *, first_attempt: bool = False) -> dict:
    if not cells or len(set(cells)) != len(cells) or any(c < 0 for c in cells):
        raise ValueError("Trace cells must be distinct zero-based IDs")
    rows = {}
    for line in content.splitlines():
        if PREFIX not in line:
            continue
        terms = line.split(PREFIX, 1)[1].split()
        data = {}
        for token in terms:
            if "=" not in token:
                raise ValueError("Malformed row trace token")
            name, value = token.split("=", 1)
            if name in data:
                raise ValueError("Duplicate field in row trace")
            data[name] = value
        if not REQUIRED.issubset(data):
            raise ValueError("Truncated row trace; missing fields")
        if data["qualified"] != "0":
            raise ValueError("Row diagnostic cannot certify a physical solution")
        it, cell = int(data["iteration"]), int(data["cell"])
        field = data["field"]
        if it <= 0 or cell not in cells or field not in {"k", "omega"}:
            raise ValueError("Invalid trace iteration, cell or field")
        if data["source_units"] != "per_unit_volume":
            raise ValueError("Source derivative units are missing or ambiguous")
        for name in ("full_residual", "closure_error", "source_identity_error", "local_relative",
                     "cell_volume", "dRk_dlogk_source", "dRk_dlogomega_source",
                     "dRk_dlogk_fv_source", "dRk_dlogomega_fv_source"):
            number = float(data[name])
            if not math.isfinite(number) or (name == "local_relative" and number < 0) or (name == "cell_volume" and number <= 0):
                raise ValueError("Non-finite or negative residual diagnostic")
        for derivative in ("dRk_dlogk", "dRk_dlogomega"):
            physical = float(data[f"{derivative}_fv_source"])
            density = float(data[f"{derivative}_source"])
            expected = density * float(data["cell_volume"])
            if not math.isclose(physical, expected, rel_tol=5e-12, abs_tol=1e-27):
                raise ValueError("Source derivative volume scaling is inconsistent")
        for name in ("production_limited", "omega_viscosity_branch", "branch_changed", "winner_accepted", "frozen_source_derivative_valid"):
            if data[name] not in ("0", "1"):
                raise ValueError("Invalid branch diagnostic flag")
        key = (it, cell, field)
        # An outer rollback can restart SST iteration numbering at 1.
        # Only an explicit --first-attempt request may segment the log.
        if (first_attempt and it == 1 and key in rows and
                any(previous[0] > 1 for previous in rows)):
            break
        if key in rows:
            raise ValueError("Duplicate row trace; ambiguous nonlinear iteration")
        rows[key] = data
    if not rows:
        raise ValueError("No V33 row traces found")
    iterations = sorted({item[0] for item in rows})
    if iterations != list(range(1, iterations[-1] + 1)):
        raise ValueError("Missing inner SST iterations")
    for it in iterations:
        for cell in cells:
            for field in ("k", "omega"):
                if (it, cell, field) not in rows:
                    raise ValueError("Missing SST field or cell in an iteration")
    summaries = {}
    for cell in cells:
        k = [rows[(i, cell, "k")] for i in iterations]
        omega = [rows[(i, cell, "omega")] for i in iterations]
        summaries[str(cell)] = {
            "k_final_residual": float(k[-1]["full_residual"]),
            "omega_final_residual": float(omega[-1]["full_residual"]),
            "k_max_relative": max(float(r["local_relative"]) for r in k),
            "omega_max_relative": max(float(r["local_relative"]) for r in omega),
            "k_final_source_dR_dlogomega_per_volume": float(k[-1]["dRk_dlogomega_source"]),
            "k_final_source_dR_dlogomega_fv": float(k[-1]["dRk_dlogomega_fv_source"]),
            "cell_volume": float(k[-1]["cell_volume"]),
            "k_source_derivative_valid_final": k[-1]["frozen_source_derivative_valid"] == "1",
            "branch_switch_iterations": [i for i, r in zip(iterations, k)
                                         if r["branch_changed"] == "1"],
        }
    return {"status": "DIAGNOSTIC_ONLY_NOT_QUALIFIED",
            "complete_inner_iterations": len(iterations),
            "cells": summaries,
            "outer_commits": "NOT_ESTABLISHED_FROM_ROW_TRACE",
            "native_openfoam_equivalence": "NOT_TESTED"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("log", type=Path)
    p.add_argument("--cells", required=True)
    p.add_argument("--first-attempt", action="store_true",
                   help="Stop at first SST iteration-number restart after a rejected outer step")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    cells = [int(v) for v in args.cells.split(",")]
    result = parse_log(args.log.read_text(errors="replace"), cells,
                       first_attempt=args.first_attempt)
    data = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
