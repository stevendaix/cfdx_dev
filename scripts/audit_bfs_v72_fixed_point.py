"""Fail-closed evidence audit for bounded PRIVATE V72 COUPLED/SST cycles.

This does not issue an industrial convergence certificate or authorize commits.
It compares the canonical V65/V69/INCOMPRESSIBLE_ITER history against V71,
then validates V72 repeated-cycle and newly reassembled six-field evidence.
"""

import argparse
import json
import math
from pathlib import Path

HISTORY = (
    "GLOBAL_SAME_OPERATOR_TRIAL_V65",
    "GLOBAL_SAME_OPERATOR_COMMITTED_V65",
    "BFS_SST_RC_PICARD_V69",
    "INCOMPRESSIBLE_ITER",
)
ATTEMPT = "BFS_JOINT_FIXED_POINT_V72"
CERT = "BFS_JOINT_CERTIFICATE_V72"
SUMMARY = "BFS_JOINT_SUMMARY_V72"


def events(text, prefix):
    return [line for line in text.splitlines() if line.startswith(prefix + " ")]


def decode(line):
    output = {}
    for word in line.split()[1:]:
        if "=" in word:
            k, v = word.split("=", 1)
            if k in output:
                raise ValueError("duplicate evidence key " + k)
            output[k] = v
    return output


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError("nonfinite " + key)
    return value


def boolean(row, key):
    if row[key] not in ("0", "1"):
        raise ValueError("invalid boolean " + key)
    return row[key] == "1"


def audit(baseline, candidate):
    counts = {}
    for prefix in HISTORY:
        reference = events(baseline, prefix)
        actual = events(candidate, prefix)
        if not reference or reference != actual:
            raise ValueError("production history differs: " + prefix)
        counts[prefix] = len(actual)
    attempts = [decode(s) for s in events(candidate, ATTEMPT)]
    certificates = [decode(s) for s in events(candidate, CERT)]
    summaries = [decode(s) for s in events(candidate, SUMMARY)]
    if not attempts or not summaries or len(summaries) != 1:
        raise ValueError("missing/duplicated V72 evidence")
    summary = summaries[0]
    if not (
        number(summary, "initial_mass") >= 0.0 and number(summary, "final_mass") >= 0.0
    ):
        raise ValueError("negative mass")
    if (
        boolean(summary, "commit_authorized")
        or boolean(summary, "committed")
        or boolean(summary, "qualified")
    ):
        raise ValueError("private V72 result unlawfully committed or qualified")
    evaluated = 0
    selected = 0
    selected_mass = []
    for row in attempts:
        if row["iteration"] != "23" or row.get("retry") != "8":
            raise ValueError("V72 probe outside specified pilot")
        for k in (
            "mass_before",
            "mass_after",
            "omega_after",
            "merit_ratio",
            "worst_growth",
            "rAU_change",
            "nut_change",
            "phi_change",
        ):
            if number(row, k) < 0.0:
                raise ValueError("negative V72 metric " + k)
        if (
            boolean(row, "committed")
            or boolean(row, "qualified")
            or not boolean(row, "private_trial")
        ):
            raise ValueError("V72 attempt published a private state")
        if boolean(row, "evaluated"):
            evaluated += 1
        if boolean(row, "selected"):
            selected += 1
            if not (
                boolean(row, "evaluated")
                and boolean(row, "v65_guards")
                and boolean(row, "mass_improved")
            ):
                raise ValueError("selected cycle violates acceptance contract")
            selected_mass.append(number(row, "mass_after"))
    if selected != int(summary["accepted_private_cycles"]) or selected != len(
        certificates
    ):
        raise ValueError("selected cycles and certificates mismatch")
    if (
        selected > int(summary["attempted_cycles"])
        or int(summary["attempted_cycles"]) > 16
    ):
        raise ValueError("invalid fixed-point budget")
    for i, row in enumerate(certificates):
        if not boolean(row, "replay"):
            raise ValueError("new F_phys replay missing")
        for k in ("face_identity", "continuity_face_gap", "rc_fixed_point"):
            if number(row, k) > 1e-10 or number(row, k) < 0.0:
                raise ValueError("nonconservative face identity or unclosed RC")
        if abs(number(row, "mass") - selected_mass[i]) > 1e-12:
            raise ValueError("certificate from a different state")
        if (
            boolean(row, "linear_momentum_reissued")
            or boolean(row, "linear_sst_reissued")
            or boolean(row, "commit_authorized")
            or boolean(row, "committed")
            or boolean(row, "qualified")
        ):
            raise ValueError("unsubstantiated production linear certificate / commit")
        if boolean(row, "mass_closed") != (number(row, "mass") <= 1e-10):
            raise ValueError("false mass closure certificate")
        for k in ("momentum_scaled", "k_scaled", "omega_scaled"):
            if number(row, k) < 0.0:
                raise ValueError("negative residual")
    if selected and abs(number(summary, "final_mass") - selected_mass[-1]) > 1e-12:
        raise ValueError("summary does not match certified final state")
    if boolean(summary, "mass_closed") != (number(summary, "final_mass") <= 1e-10):
        raise ValueError("false summary mass closure")
    return {
        "status": "diagnostic_only_no_commit_certificate",
        "historical_identity": counts,
        "attempts": len(attempts),
        "evaluated": evaluated,
        "accepted_private_cycles": selected,
        "final_physical_mass": number(summary, "final_mass"),
        "mass_closed": boolean(summary, "mass_closed"),
        "full_production_convergence_qualified": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.baseline.read_text(), args.log.read_text())
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
