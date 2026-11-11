#!/usr/bin/env python3
"""Fail-closed first-attempt V26 log-Jv/L2 audit, never a BFS certification."""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path


def kv(line: str) -> dict[str, str]:
    return dict(re.findall(r"([A-Za-z_][A-Za-z_0-9]*)=([^\s]+)", line))


def first_attempt(log: str) -> dict:
    first = re.split(r"CFDX nonlinear iteration 1 rejected|retry=1", log, maxsplit=1)[0]
    rows = [kv(row) for row in first.splitlines() if "SST_JFNK_V22 iteration=" in row]
    candidates = [row for row in rows if row.get("v26_active") == "1" and
                  row.get("v26_jvp_available") == "1"]
    if not candidates:
        raise ValueError("no V26 multiscale audit on first physical attempt")
    row = candidates[-1]
    def number(key: str) -> float:
        try:
            val = float(row[key])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"missing V26 key: {key}") from exc
        if not math.isfinite(val):
            raise ValueError(f"nonfinite V26 key: {key}")
        return val
    steps = [number(f"v26_jvp_h{i}") for i in (1,2,3)]
    if not (steps[0] > steps[1] > steps[2] > 0):
        raise ValueError("non-decreasing difference scale")
    mismatch = number("v26_jvp_max_disagreement")
    if mismatch < 0 or number("v26_jvp_evaluations") != 7:
        raise ValueError("invalid log Jacobian audit")
    baseline_l2 = number("v26_l2_before")
    best_trial = number("v26_l2_best_trial")
    baseline_linf = number("old_merit")
    selected_linf = number("new_merit")
    if baseline_l2 <= 0 or baseline_linf <= 0 or best_trial < 0 or selected_linf < 0:
        raise ValueError("invalid nonlinear merit")
    smooth = row.get("v26_l2_accepted") == "1"
    selected = row.get("selected")
    if selected not in {"picard", "jfnk", "jfnk_block", "jfnk_spatial",
                        "jfnk_spatial_log_v25", "jfnk_smooth_v26", "none"}:
        raise ValueError("unknown selected solver branch")
    if smooth:
        alpha = number("v26_alpha")
        l2 = number("v26_l2_selected")
        physical = number("v26_linf_selected")
        if not (0 < alpha <= 1 and l2 < baseline_l2*(1-1e-9*alpha) and
                physical < baseline_linf*(1-1e-9*alpha)):
            raise ValueError("V26 would promote L2-only or non-descent step")
    elif selected == "jfnk_smooth_v26":
        raise ValueError("selected unavailable smooth candidate")
    if row.get("qualified") == "1":
        raise ValueError("diagnostic row must never assert BFS certification")
    if selected != "none" and not selected_linf < baseline_linf:
        raise ValueError("final selection did not reduce true physical Linf")
    return {
        "jvp_perturbations": steps,
        "maximum_directional_difference": mismatch,
        "l2_before": baseline_l2,
        "l2_best_trial": best_trial,
        "l2_proposal_physically_admissible": smooth,
        "selected_solver": selected,
        "physical_linf_before": baseline_linf,
        "physical_linf_after": selected_linf,
        "sst_converged": False,
        "bfs_converged": False,
        "nasa_validated": False,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log", type=Path)
    args = ap.parse_args()
    print(json.dumps(first_attempt(args.log.read_text()), indent=2))


if __name__ == "__main__":
    main()
