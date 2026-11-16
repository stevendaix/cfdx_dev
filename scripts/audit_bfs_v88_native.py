"""Fail-closed audit of REAL V88 pre-SST local-step decisions and V85 commits."""

import argparse
import json
import math
import re
from pathlib import Path

MARKER = re.compile(r"\b([a-zA-Z_]+)=([^\s]+)")


def fields(line):
    return dict(MARKER.findall(line))


def audit(log_text):
    trials, committed, failures = [], [], []
    for lineno, line in enumerate(log_text.splitlines(), 1):
        if line.startswith("BFS_LOCAL_STEP_TRIAL_V88 "):
            f = fields(line)
            required = (
                "iteration",
                "retry",
                "continue_sst",
                "max_delta_u_over_inlet",
                "delta_u_limit",
                "worst_cell",
                "delta_pressure_span_over_qin",
                "alpha_u",
                "alpha_p",
                "reason",
                "committed",
                "qualified",
            )
            if any(key not in f for key in required):
                raise ValueError(f"Incomplete V88 trial at line {lineno}")
            if f["committed"] != "0" or f["qualified"] != "0":
                raise ValueError("Pre-SST trial falsely certified")
            r = {
                "line": lineno,
                "iteration": int(f["iteration"]),
                "retry": int(f["retry"]),
                "continue_sst": int(f["continue_sst"]),
                "max_delta_u_over_inlet": float(f["max_delta_u_over_inlet"]),
                "delta_u_limit": float(f["delta_u_limit"]),
                "worst_cell": int(f["worst_cell"]),
                "delta_pressure_span_over_qin": float(
                    f["delta_pressure_span_over_qin"]
                ),
                "alpha_u": float(f["alpha_u"]),
                "alpha_p": float(f["alpha_p"]),
                "reason": f["reason"],
            }
            if not all(
                math.isfinite(r[key]) and r[key] >= 0
                for key in (
                    "max_delta_u_over_inlet",
                    "delta_u_limit",
                    "delta_pressure_span_over_qin",
                    "alpha_u",
                    "alpha_p",
                )
            ):
                raise ValueError("Nonfinite V88 evidence")
            expected = int(r["max_delta_u_over_inlet"] <= r["delta_u_limit"] + 1e-10)
            if r["continue_sst"] not in (0, 1) or r["continue_sst"] != expected:
                raise ValueError("V88 verdict inconsistent with recorded increment")
            if (expected == 0) != (r["reason"] == "local_velocity_increment_envelope"):
                raise ValueError("V88 reason inconsistent")
            trials.append(r)
        elif line.startswith("GLOBAL_PROGRESS_COMMITTED_V85 "):
            f = fields(line)
            if f.get("qualified") != "0" and f.get("converged") != "1":
                raise ValueError("Unqualified provisional commit mislabeled")
            if "iteration" in f:
                committed.append(int(f["iteration"]))
        elif line.startswith("CFDX nonlinear iteration ") and "rejected;" in line:
            failures.append((lineno, line[:300]))
    if committed != list(range(1, len(committed) + 1)):
        raise ValueError("Committed iteration sequence non-contiguous")
    # A V88-rejected trial must not enter SST before the retry of that iteration.
    lines = log_text.splitlines()
    for r in trials:
        if r["continue_sst"] == 0:
            for line in lines[r["line"] :]:
                if line.startswith(
                    (
                        "BFS_LOCAL_STEP_TRIAL_V88 ",
                        "BFS_EARLY_FIELD_TRIAL_V87 ",
                        "CFDX nonlinear iteration ",
                    )
                ):
                    break
                if "phase=SST_TRANSPORT_BEGIN" in line:
                    raise ValueError("V88 rejected trial nevertheless entered SST")
    return {
        "schema": "CFDX_V88_NATIVE_AUDIT_1",
        "trial_count": len(trials),
        "pre_sst_rejected": sum(1 - r["continue_sst"] for r in trials),
        "accepted_outer_count": len(committed),
        "accepted_outer_indices": committed,
        "trials": trials,
        "retries_in_log": len(failures),
        "finally_converged": False,
        "physical_qualification": "NOT_QUALIFIED",
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("log", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    report = audit(a.log.read_text(errors="replace"))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "trials"}, indent=2))
