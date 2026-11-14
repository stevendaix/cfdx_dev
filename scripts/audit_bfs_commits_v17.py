"""Fail-closed accepted-state ledger for BFS U,p,k,omega.

A V13 committed record is useful diagnostic evidence, but *does not* certify
same-operator contraction: only the V14 ledger carries that contract. Neither
ledger certifies a physically correct solution or agreement with NASA.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

try:
    from qualify_bfs_global_contraction_v14 import certify, QualificationError
except ImportError:
    from scripts.qualify_bfs_global_contraction_v14 import certify, QualificationError

SCHEMA = "CFDX-BFS-COMMITTED-STATES/v17"


def _parse(line: str) -> dict:
    values = dict(re.findall(r"(\w+)=([^\s]+)", line))
    mandatory = ("iteration", "U", "p_continuity", "k", "omega", "nut_closure",
                 "merit", "converged", "qualified")
    if any(key not in values for key in mandatory):
        raise ValueError("incomplete accepted-state certificate")
    try:
        iteration = int(values["iteration"])
        residuals = {name: float(values[name]) for name in mandatory[1:7]}
        converged = int(values["converged"])
        qualified = int(values["qualified"])
    except (ValueError, OverflowError) as exc:
        raise ValueError("malformed accepted-state certificate") from exc
    if iteration <= 0 or converged not in (0, 1) or qualified != 0 or any(
        v < 0 or not math.isfinite(v) for v in residuals.values()
    ):
        raise ValueError("nonfinite, negative or unsupported qualified state")
    return dict(iteration=iteration, **residuals, converged=bool(converged))


def audit(text: str) -> dict:
    states = [_parse(line) for line in text.splitlines()
              if line.startswith("GLOBAL_NONLINEAR_COMMITTED_V13 ")]
    if any(next_state["iteration"] <= state["iteration"]
           for state, next_state in zip(states, states[1:])):
        raise ValueError("out-of-order or duplicate committed iterations")
    record = {"schema": SCHEMA, "committed_states": len(states),
              "last_committed": states[-1] if states else None,
              "three_states_observed": len(states) >= 3,
              "strict_same_operator_contraction": False,
              "fully_converged": False,
              "physical_qualified": False,
              "reference_comparison_permitted": False,
              "status": "NOT_QUALIFIED", "states": states}
    try:
        v14 = certify(text, minimum_commits=3)
    except QualificationError as exc:
        record["contraction_reason"] = str(exc)
    else:
        # An opt-in V14 ledger, not simply decreasing V13 merits, is necessary
        # to demonstrate three accepted contractions with one operator.
        record["strict_same_operator_contraction"] = bool(v14["strict_contraction_proved"])
        record["fully_converged"] = bool(v14["physical_converged"])
        record["contraction_certificate"] = v14
        # Even a numerically converged V14 ledger does NOT certify physical
        # accuracy without independent NASA, wall, and mesh studies.
    if not states:
        record["first_blocker"] = "NO_COMMITTED_U_P_K_OMEGA_STATE"
    elif len(states) < 3:
        record["first_blocker"] = "FEWER_THAN_THREE_COMMITTED_STATES"
    elif not record["strict_same_operator_contraction"]:
        record["first_blocker"] = "SAME_OPERATOR_CONTRACTION_LEDGER_MISSING_OR_FAILED"
    elif not record["fully_converged"]:
        record["first_blocker"] = "NONLINEAR_PHYSICAL_RESIDUAL_NOT_CLOSED"
    else:
        record["first_blocker"] = "NASA_WALL_MESH_QUALIFICATION_MISSING"
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(args.log.read_text(encoding="utf-8", errors="replace"))
        code = 0 if result["strict_same_operator_contraction"] and result["fully_converged"] else 3
    except (OSError, ValueError) as exc:
        result = {"schema": SCHEMA, "status": "INVALID_EVIDENCE",
                  "first_blocker": str(exc), "physical_qualified": False,
                  "reference_comparison_permitted": False}
        code = 3
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
