#!/usr/bin/env python3
"""Strict V70 BFS152 native-log replay and fail-closed staging audit.

This tool checks evidence, not CFD convergence. Original V69 accepted and
rejected transactions MUST be byte-identical when V70's private probes run.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

HISTORICAL_TAGS = (
    "GLOBAL_SAME_OPERATOR_TRIAL_V65",
    "GLOBAL_SAME_OPERATOR_COMMITTED_V65",
    "BFS_SST_RC_PICARD_V69",
    "INCOMPRESSIBLE_ITER",
)


def records(text: str, tag: str) -> list[str]:
    return [line for line in text.splitlines() if line.startswith(tag + " ")]


def parse_values(line: str) -> dict[str, str]:
    values = {}
    for entry in line.split()[1:]:
        if "=" in entry:
            key, value = entry.split("=", 1)
            values[key] = value
    return values


def numeric(values: dict[str, str], name: str) -> float:
    number = float(values[name])
    if not math.isfinite(number):
        raise ValueError(f"nonfinite {name}")
    return number


def check_v70_evidence(log: str, reference: str) -> dict:
    counts = {}
    for tag in HISTORICAL_TAGS:
        original = records(reference, tag)
        current = records(log, tag)
        if not original or original != current:
            raise ValueError(f"{tag}: missing reference or changed production trajectory")
        counts[tag] = len(current)
    if counts["GLOBAL_SAME_OPERATOR_TRIAL_V65"] != 149 or \
       counts["GLOBAL_SAME_OPERATOR_COMMITTED_V65"] != 22 or \
       counts["BFS_SST_RC_PICARD_V69"] != 6:
        raise ValueError("incomplete historical BFS152/V69 replay")
    entries = records(log, "BFS_SST_RC_CONSERVATIVE_V70")
    if not entries:
        raise ValueError("V70 pressure/flux projection evidence absent")
    closed_count = 0
    trials = []
    for row in entries:
        v = parse_values(row)
        if v.get("committed") != "0" or v.get("qualified") != "0":
            raise ValueError("private V70 trial was published or falsely qualified")
        if v.get("iteration") != "23" or v.get("retry") != "8":
            raise ValueError("V70 tested a different outer state")
        if "cycle" in v:
            for key in ("fraction", "mass_before", "mass_frozen_after",
                        "mass_true_after_sst", "pressure_delta", "face_delta",
                        "merit_ratio", "worst_growth"):
                numeric(v, key)
            frac = numeric(v, "fraction")
            frozen_before = numeric(v, "mass_before")
            frozen_after = numeric(v, "mass_frozen_after")
            true_after = numeric(v, "mass_true_after_sst")
            if not (0 < frac <= 1):
                raise ValueError("V70 invalid pressure continuation fraction")
            # A fractional projection is intentionally NOT a certified
            # conservative correction. It must only decrease the frozen
            # mass imbalance, and must not be mislabelled as a closed step.
            if frac == 1 and frozen_after > 1e-10:
                raise ValueError("V70 full projection did not conserve frozen mass")
            if frac < 1 and not frozen_after < frozen_before:
                raise ValueError("V70 partial projection did not reduce frozen divergence")
            if v.get("candidate_eligible") == "1" and true_after > 1e-10:
                raise ValueError("V70 false final mass-certification claim")
            closed_count += 1
            trials.append({k: v[k] for k in (
                "cycle", "fraction", "mass_before", "mass_frozen_after", "mass_true_after_sst",
                "pressure_delta", "face_delta", "merit_ratio", "worst_growth",
                "candidate_eligible")})
        elif v.get("closed") != "0" and v.get("certified_private") != "1":
            raise ValueError("V70 unexplained diagnostic record")
    tx = records(log, "BFS_SST_RC_TRANSACTION_V70")
    if len(tx) != 1:
        raise ValueError("V70 rollback/nonmutation proof absent")
    values = parse_values(tx[0])
    if (values.get("U_p_phi_k_omega_nut_unchanged") != "1" or
        values.get("committed") != "0" or values.get("qualified") != "0"):
        raise ValueError("V70 carrier/SST snapshot changed")
    return {
        "status": "forensic_only_not_qualified",
        "historical_counts": counts,
        "v70_projection_cycles": closed_count,
        "v70_rejection_rows": len(entries) - closed_count,
        "v70_projection_rows": trials,
        "outer_transaction_changed": False,
        "real_v70_commit": False,
        "physical_bfs_converged": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = check_v70_evidence(args.log.read_text(), args.baseline.read_text())
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())