#!/usr/bin/env python3
"""Classify the *first observed* BFS blocker using production stdout only.

A timeout is an incomplete observation, never an assertion that a solver
algorithm failed or that the experimental physics did not converge.  The
reported physical stage comes from opt-in BFS_PHASE_V10 markers emitted by the
C++ production implementation, not the requested Schur model name.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

_PHASE = re.compile(r"\bBFS_PHASE_V10\s+iteration=(\d+).*?\bphase=([A-Z_]+)")
_ACTUAL_LINEAR_PLAN = re.compile(
    r"\bBFS_PHASE_V10\s+iteration=(\d+)\s+phase=COUPLED_LINEAR_END.*?"
    r"\bactual_krylov=([a-z0-9_]+)\s+actual_preconditioner=([a-z0-9_]+)"
    r"\s+resolved_schur=([a-z0-9_]+)")
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_LINEAR_FAILURE = re.compile(
    r"(pressure-correction|coupled momentum-continuity) solve did not converge"
    r" \(status=(\d+), iterations=(\d+), residual=(" + _NUMBER +
    r"), relative=(" + _NUMBER + r")\)")
_SST_FAILURE = re.compile(r"SST update rejected: ([^\r\n]+)")
_GLOBAL_FAILURE = re.compile(r"(global nonlinear|post-SST physical momentum|SST residual/boundedness gate failed)[^\r\n]*", re.I)


def classify_log(content: str, *, execution: str, exit_code: int | None) -> dict:
    """Fail closed; a stopped run has only an observed stage, not convergence."""
    phases = [(m.group(2), int(m.group(1))) for m in _PHASE.finditer(content)]
    resolved_plans = [
        {"iteration": int(m.group(1)), "actual_krylov": m.group(2),
         "actual_preconditioner": m.group(3), "resolved_schur": m.group(4)}
        for m in _ACTUAL_LINEAR_PLAN.finditer(content)
    ]
    events: list[tuple[int, str, str, dict]] = []
    for match in _LINEAR_FAILURE.finditer(content):
        stage = "PRESSURE" if match.group(1) == "pressure-correction" else "KRYLOV_SCHUR"
        events.append((match.start(), stage, match.group(0), {
            "solver_status_code": int(match.group(2)),
            "linear_iterations": int(match.group(3)),
            "true_residual": float(match.group(4)),
            "true_residual_relative": float(match.group(5)),
        }))
    for match in _SST_FAILURE.finditer(content):
        events.append((match.start(), "SST", match.group(0), {}))
    for match in _GLOBAL_FAILURE.finditer(content):
        events.append((match.start(), "GLOBAL_U_P_K_OMEGA", match.group(0), {}))
    events.sort(key=lambda e: e[0])
    result: dict = {
        "schema": "CFDX-BFS-FIRST-BLOCKER/v10",
        "first_blocker": "UNRESOLVED",
        "first_failure_message": None,
        "last_observed_phase": phases[-1][0] if phases else None,
        "phase_marker_count": len(phases),
        "first_resolved_linear_plan": resolved_plans[0] if resolved_plans else None,
        "resolved_linear_plan_observations": len(resolved_plans),
        "observed_nonlinear_iteration": phases[-1][1] if phases else None,
        "linear_evidence": None,
        "proof_status": "NOT_QUALIFIED",
    }
    if events:
        _, result["first_blocker"], result["first_failure_message"], details = events[0]
        result["linear_evidence"] = details or None
        return result
    if execution == "TIMEOUT_INCOMPLETE":
        last = result["last_observed_phase"]
        stage = {
            "PRESSURE_LINEAR_BEGIN": "PRESSURE_RUNNING_TIMEOUT",
            "COUPLED_LINEAR_BEGIN": "KRYLOV_SCHUR_RUNNING_TIMEOUT",
            "SST_TRANSPORT_BEGIN": "SST_RUNNING_TIMEOUT",
        }.get(last, "UNRESOLVED_SETUP_OR_ASSEMBLY_TIMEOUT")
        result["first_blocker"] = stage
        return result
    if exit_code == 0:
        result["first_blocker"] = "NO_FAILURE_REPORTED_REQUIRES_INDEPENDENT_GATE"
    elif "INCOMPRESSIBLE_ITER" in content:
        result["first_blocker"] = "POST_MOMENTUM_NOT_CERTIFIED"
    else:
        result["first_blocker"] = "UNRESOLVED_START_OR_ASSEMBLY"
    return result


def diagnose_campaign(path: Path) -> dict:
    campaign = json.loads(path.read_text(encoding="utf-8"))
    cases = []
    for item in campaign.get("cases", []):
        if not str(item.get("case", "")).startswith("bfs_"):
            continue
        log_file = Path(item["log"])
        if not log_file.is_file():
            diagnostic = {"first_blocker": "LOG_MISSING", "proof_status": "NOT_QUALIFIED"}
            sha = None
        else:
            contents = log_file.read_bytes()
            sha = hashlib.sha256(contents).hexdigest()
            diagnostic = classify_log(contents.decode("utf-8", "replace"),
                execution=item.get("execution", "NOT_RUN"), exit_code=item.get("exit_code"))
        cases.append({
            "case": item["case"],
            "execution": item.get("execution"),
            "exit_code": item.get("exit_code"),
            "wall_seconds": item.get("wall_seconds"),
            "case_hdf5_sha256": item.get("mesh_sha256"),
            "binary_sha256": item.get("binary_sha256"),
            "log_sha256": sha,
            "nonlinear_v8_proven": item.get("bfs_v8_gate", {}).get("nonlinear_v8_proven", False),
            "independent_observation": diagnostic,
        })
    return {"schema": "CFDX-BFS-DIAGNOSTIC-MATRIX/v10",
            "physical_qualification": "NOT_EVALUATED",
            "cases": cases}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = diagnose_campaign(args.campaign)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for item in record["cases"]:
        print(item["case"], item["execution"],
              item["independent_observation"]["first_blocker"])
    return 0 if record["cases"] and all(
        x["nonlinear_v8_proven"] for x in record["cases"]) else 2


if __name__ == "__main__":
    raise SystemExit(main())
