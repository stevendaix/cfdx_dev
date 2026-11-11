#!/usr/bin/env python3
"""Fail-closed V27 first-attempt diagnostic. Never certifies BFS/NASA physics."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from diagnose_bfs_sst_v26 import kv


def first_attempt(log: str) -> dict:
    # The first attempted outer state is the only comparable one. A retry
    # may use a modified damping configuration and is NOT the same experiment.
    first = log.split("CFDX nonlinear iteration 1 rejected", 1)[0]
    traces = [kv(line) for line in first.splitlines()
              if "SST_JFNK_V22 iteration=10" in line]
    if not traces:
        raise ValueError("no V27 first-attempt iteration-10 trace")
    trace = traces[0]
    if trace.get("v27_active") != "1" or trace.get("v26_active") != "1":
        raise ValueError("V27 opt-in not present in tested first attempt")
    if trace.get("qualified") != "0":
        raise ValueError("diagnostic cannot qualify a BFS solution")
    def number(field: str) -> float:
        try:
            n = float(trace[field])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"V27 required metric missing: {field}") from exc
        if not math.isfinite(n):
            raise ValueError(f"V27 nonfinite required metric: {field}")
        return n
    before = number("old_merit")
    after = number("new_merit")
    active = trace.get("v27_accepted") == "1"
    selected = trace.get("selected", "")
    if selected not in ("none", "picard", "jfnk", "jfnk_block",
                        "jfnk_spatial", "jfnk_spatial_log_v25",
                        "jfnk_smooth_v26", "sst_subspace_trust_v27"):
        raise ValueError("unknown selected branch")
    if selected == "sst_subspace_trust_v27" and not active:
        raise ValueError("unavailable V27 step selected")
    if selected != "none" and not (after < before):
        raise ValueError("selected step has no true physical Linf descent")
    attempts = number("v27_evaluations")
    bases = number("v27_bases")
    if not (bases in (0, 3, 4) and attempts >= 0):
        raise ValueError("V27 basis/evaluation accounting invalid")
    if active:
        proposed = number("v27_linf")
        l2 = number("v27_l2")
        old_l2 = number("v26_l2_before")
        radius = number("v27_radius")
        ratio = number("v27_ratio")
        slope = number("v27_slope")
        if not (0 < radius <= 1.5 and ratio > 1e-3 and slope < 0 and
                l2 < old_l2*(1-1e-9*radius) and
                proposed < before*(1-1e-9*radius)):
            raise ValueError("V27 candidate fails one of its physical gates")
        if trace.get("v27_choice") not in {
            "projected_l2_cauchy", "active_linf_restoration"}:
            raise ValueError("unsupported V27 direction type")
    return {
        "v27_accepted": active,
        "v27_selected": selected == "sst_subspace_trust_v27",
        "v27_evaluations": int(attempts),
        "v27_bases": int(bases),
        "physical_linf_before": before,
        "physical_linf_after": after,
        "sst_converged": False,
        "bfs_converged": False,
        "nasa_validated": False,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("log",type=Path)
    args=ap.parse_args()
    print(json.dumps(first_attempt(args.log.read_text()), indent=2))

if __name__ == "__main__":
    main()
