"""Fail-closed BFS152 V65 evaluation of same-state V47 versus historical COUPLED.

This is an execution/cause audit, NOT physical validation. A solver exit with
'Converged NO' is never promoted to PASS by a falling merit or a CTest PASS.
The physical model must be identical except explicitly listed globalization
flags, and every successful V65 transaction must be traceable to a same-state
Fphys trial (old and candidate RC/rAU, six blocks, common frozen scales).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

MARKER = "GLOBAL_SAME_OPERATOR_TRIAL_V65 "
COMMIT = "GLOBAL_SAME_OPERATOR_COMMITTED_V65 "
GLOBAL = "GLOBAL_NONLINEAR_COMMITTED_V13 "
KVS = re.compile(r"([A-Za-z][A-Za-z_0-9]*)=([^\s]+)")


def fields(line):
    return dict(KVS.findall(line))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(text, prefix):
    out = []
    for line in text.splitlines():
        if line.startswith(prefix):
            out.append(fields(line))
    return out


def checked_number(row, key, nonnegative=True):
    if key not in row:
        raise ValueError(f"missing {key}")
    v = float(row[key])
    if not math.isfinite(v) or (nonnegative and v < 0):
        raise ValueError(f"nonfinite or negative {key}")
    return v


def checked_int(row, key, choices=None):
    if key not in row:
        raise ValueError(f"missing {key}")
    v = int(row[key])
    if choices is not None and v not in choices:
        raise ValueError(f"bad {key}")
    return v


def physical_rows(log):
    result = {}
    for d in rows(log, GLOBAL):
        i = checked_int(d, "iteration")
        if i < 1 or i in result:
            raise ValueError("duplicate/invalid physical commit")
        for k in ("U", "p_continuity", "k", "omega", "merit"):
            checked_number(d, k)
        checked_int(d, "qualified", {0})
        checked_int(d, "converged", {0, 1})
        result[i] = {
            k: float(d[k]) for k in ("U", "p_continuity", "k", "omega", "merit")
        }
        result[i]["converged"] = bool(int(d["converged"]))
    if not result or sorted(result) != list(range(1, len(result) + 1)):
        raise ValueError("missing/fragmented physical committed states")
    return result


def old_baseline_at_iteration(log, iteration):
    if log is None:
        return None
    physical = physical_rows(log)
    if iteration not in physical:
        return None
    return {
        "iteration": iteration,
        "U": physical[iteration]["U"],
        "continuity": physical[iteration]["p_continuity"],
    }


def audit(log, baseline_log=None, historical_v47_log=None):
    trials = {}
    for d in rows(log, MARKER):
        i = checked_int(d, "iteration")
        retry = checked_int(d, "retry")
        if i < 1 or retry < 0 or (i, retry) in trials:
            raise ValueError("duplicate/bad trial")
        checked_number(d, "beta")
        for k in (
            "baseline",
            "joint_ok",
            "monotone_ok",
            "accepted",
            "same_state_operator",
            "qualified",
            "committed",
        ):
            checked_int(d, k, {0, 1})
        if (
            d["same_state_operator"] != "1"
            or d["qualified"] != "0"
            or d["committed"] != "0"
        ):
            raise ValueError("fabricated same-state/physical qualification")
        for k in (
            "old_merit",
            "trial_merit",
            "ratio",
            "largest_growth",
            "old_rc_gap",
            "trial_rc_gap",
            "old_div_gap",
            "trial_div_gap",
            "old_authoritative_phi_gap",
            "trial_authoritative_phi_gap",
        ):
            checked_number(d, k)
        if (
            checked_number(d, "old_rc_gap") > 1e-10
            or checked_number(d, "trial_rc_gap") > 1e-10
            or checked_number(d, "old_div_gap") > 1e-10
            or checked_number(d, "trial_div_gap") > 1e-10
        ):
            raise ValueError("unclosed independent flux / continuity")
        if d["baseline"] == "1" and i != 1:
            raise ValueError("illicit later baseline")
        if (
            d["accepted"] == "1"
            and i > 1
            and (
                d["joint_ok"] != "1"
                or d["monotone_ok"] != "1"
                or float(d["ratio"]) > 0.995 + 1e-12
            )
        ):
            raise ValueError("V65 accepted a trial that failed the existing guard")
        trials[i, retry] = d
    commits = {}
    for d in rows(log, COMMIT):
        i = checked_int(d, "iteration")
        if i < 1 or i in commits:
            raise ValueError("duplicated V65 commit")
        checked_int(d, "baseline", {0, 1})
        checked_int(d, "physical_converged", {0, 1})
        checked_int(d, "qualified", {0})
        checked_number(d, "merit")
        checked_number(d, "ratio")
        commits[i] = d
    if not commits or sorted(commits) != list(range(1, len(commits) + 1)):
        raise ValueError("no consecutive V65 committed history")
    physical = physical_rows(log)
    if set(commits) != set(physical):
        raise ValueError("missing physical certificate/commit")
    records = []
    for i, commit in sorted(commits.items()):
        accepted = [
            t for (it, _), t in trials.items() if it == i and t["accepted"] == "1"
        ]
        if len(accepted) != 1:
            raise ValueError("missing/duplicated accepted trial for commit")
        trial = accepted[0]
        for key in ("ratio", "trial_merit"):
            mark = "merit" if key == "trial_merit" else key
            if not math.isclose(
                float(trial[key]), float(commit[mark]), rel_tol=1e-10, abs_tol=1e-13
            ):
                raise ValueError("trial/commit mismatch")
        if bool(int(commit["physical_converged"])) != physical[i]["converged"]:
            raise ValueError("physical qualification not consistent")
        records.append(
            {
                "iteration": i,
                "beta": float(trial["beta"]),
                "retry": int(trial["retry"]),
                "same_operator_merit": float(trial["trial_merit"]),
                "merit_ratio": float(trial["ratio"]),
                "U_post_sst": physical[i]["U"],
                "continuity": physical[i]["p_continuity"],
            }
        )
    if any(
        int(t["accepted"]) and int(t["iteration"]) not in commits
        for t in trials.values()
    ):
        raise ValueError("uncommitted accepted V65 candidate")
    base = None
    historical = None
    if baseline_log is not None:
        old = physical_rows(baseline_log)
        base = {
            "committed_steps": len(old),
            "last_momentum": old[max(old)]["U"],
            "last_continuity": old[max(old)]["p_continuity"],
            "worsening_transitions": sum(
                old[i]["U"] > old[i - 1]["U"] for i in range(2, len(old) + 1)
            ),
        }
    if historical_v47_log is not None:
        old = physical_rows(historical_v47_log)
        historical = {
            "committed_steps": len(old),
            "last_momentum": old[max(old)]["U"],
            "last_continuity": old[max(old)]["p_continuity"],
            "failed_at_iteration": len(old) + 1,
        }
    # A printed terminal marker without 40 actual transaction certificates is
    # insufficient: an incomplete, forged, or truncated log must fail closed.
    complete = len(records) == 40 and (
        "Converged NO Iterations 40" in log or "Converged YES Iterations 40" in log
    )
    failed_index = len(records) + 1
    failed_trials = [d for (i, _), d in sorted(trials.items()) if i == failed_index]
    halted = {
        "iteration": failed_index,
        "trials": len(failed_trials),
        "all_rejected": bool(failed_trials)
        and all(d["accepted"] == "0" for d in failed_trials),
        "smallest_beta": min((float(d["beta"]) for d in failed_trials), default=None),
        "last_ratio": float(failed_trials[-1]["ratio"]) if failed_trials else None,
        "last_joint_ok": int(failed_trials[-1]["joint_ok"]) if failed_trials else None,
        "last_monotone_ok": int(failed_trials[-1]["monotone_ok"])
        if failed_trials
        else None,
    }
    has_final_convergence = "Converged YES" in log and any(
        x["converged"] for x in physical.values()
    )
    return {
        "schema": "CFDX-BFS-V65-COMPARABLE-MERIT/1",
        "comparison": "same state-only Fphys at old and candidate; old-state block scaling frozen",
        "iterations_committed": len(records),
        "trial_count": len(trials),
        "rejected_trials": sum(d["accepted"] == "0" for d in trials.values()),
        "min_beta": min(x["beta"] for x in records),
        "last_U": records[-1]["U_post_sst"],
        "last_continuity": records[-1]["continuity"],
        "last_merit": records[-1]["same_operator_merit"],
        "complete_40_iteration_run": complete,
        "certified_convergence": has_final_convergence,
        "physical_qualification": False,
        "historical_coupled": base,
        "historical_v47": historical,
        "halted_trial": halted,
        "historical_at_same_iteration": (
            old_baseline_at_iteration(baseline_log, len(records)) if base else None
        ),
        "max_authoritative_flux_disagreement": max(
            max(
                float(t["old_authoritative_phi_gap"]),
                float(t["trial_authoritative_phi_gap"]),
            )
            for t in trials.values()
        ),
        "records": records,
    }


def physical_setup_equal(baseline, corrected):
    import h5py

    with h5py.File(baseline) as h, h5py.File(corrected) as g:
        if (
            h.attrs["mesh_hash"] != g.attrs["mesh_hash"]
            or h.attrs["topology_hash"] != g.attrs["topology_hash"]
            or h.attrs["geometry_hash"] != g.attrs["geometry_hash"]
        ):
            raise ValueError("BFS mesh/geometry mismatch")
        a = json.loads(h.attrs["case_setup_json"])
        b = json.loads(g.attrs["case_setup_json"])
        flags = (
            "v13_joint_globalization",
            "v14_monotone_global_merit",
            "coupled_composite_backtracking_v47_experimental",
            "coupled_composite_backtracking_retries_v47",
        )
        for k in flags:
            b["numerics"]["raw_settings"].pop(k, None)
        if a != b:
            raise ValueError("physical case changed beyond experimental globalization")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--v65-log", required=True, type=Path)
    p.add_argument("--baseline-log", required=True, type=Path)
    p.add_argument("--historical-v47-log", required=True, type=Path)
    p.add_argument("--baseline-case", required=True, type=Path)
    p.add_argument("--v65-case", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--allow-partial", action="store_true")
    args = p.parse_args()
    physical_setup_equal(args.baseline_case, args.v65_case)
    result = audit(
        args.v65_log.read_text(),
        args.baseline_log.read_text(),
        args.historical_v47_log.read_text(),
    )
    result["sha256"] = {
        k: sha(p)
        for k, p in [
            ("v65_log", args.v65_log),
            ("baseline_log", args.baseline_log),
            ("baseline_case", args.baseline_case),
            ("v65_case", args.v65_case),
        ]
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "V65_SAME_OPERATOR_AUDIT steps=",
        result["iterations_committed"],
        " trials=",
        result["trial_count"],
        " complete=",
        int(result["complete_40_iteration_run"]),
        " physically_qualified=0",
    )
    if not result["complete_40_iteration_run"] and not args.allow_partial:
        raise ValueError(
            "V65 40-iteration execution was incomplete; report is diagnostic only"
        )


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, RuntimeError) as ex:
        raise SystemExit("V65_REJECTED: " + str(ex))
