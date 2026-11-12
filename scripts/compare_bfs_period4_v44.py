"""Compare two real production logs at committed-state level (V44 vs V43).

A four-phase predictor is *not* a new physical oracle. We require matching
frozen initial COUPLED A,b,x0 and identical committed merits before the first
full-state intervention. After that intervention, report both trajectories
without labelling either as physically qualified. No raw Krylov residual can
substitute for a committed reassembled U/p/k/omega defect.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

_C = re.compile(r"^GLOBAL_NONLINEAR_COMMITTED_V13\s+(.*)$")
_F = re.compile(r"^COUPLED_MATRIX_FINGERPRINT_V40\s+(.*)$")
_P = re.compile(r"^COUPLED_PERIOD4_FULL_V44\s+(.*)$")
_S = re.compile(r"^OUTER_FIELD_PERIOD_V43\s+(.*)$")
_Q = re.compile(r"^Converged (YES|NO) Iterations (\d+)$")


def parse_pairs(value: str) -> dict[str, str]:
    values = {}
    for token in value.split():
        if "=" in token:
            name, item = token.split("=", 1)
            if name in values:
                raise ValueError(f"duplicate {name}")
            values[name] = item
    return values


def extract(text: str) -> dict:
    commits: dict[int, float] = {}
    states: dict[int, dict] = {}
    initial = None
    triggers = []
    final = None
    for line in text.splitlines():
        if m := _F.match(line):
            if initial is not None:
                raise ValueError("multiple matrix fingerprints")
            values = parse_pairs(m.group(1))
            initial = {k: values[k] for k in ("A", "b", "x0", "rows", "nnz")}
        if m := _C.match(line):
            values = parse_pairs(m.group(1))
            n = int(values["iteration"])
            merit = float(values["merit"])
            if n in commits or not math.isfinite(merit) or merit < 0:
                raise ValueError("duplicated / invalid committed merit")
            if values.get("qualified") != "0":
                raise ValueError("unjustified global physical certification")
            if values.get("converged") != "0":
                raise ValueError("converged state needs separate complete qualification")
            if any(not math.isfinite(float(values[k])) or float(values[k]) < 0
                   for k in ("U", "p_continuity", "k", "omega", "nut_closure")):
                raise ValueError("invalid physical component residual")
            commits[n] = merit
        if m := _P.match(line):
            values = parse_pairs(m.group(1))
            n = int(values["iteration"])
            if (values.get("k_omega_averaged") != "1" or
                values.get("nut_recomputed") != "1" or
                values.get("strain_recomputed") != "1" or
                values.get("committed") != "0" or
                values.get("physical_qualified") != "0"):
                raise ValueError("invalid full-phase predictor disclosure")
            triggers.append(n)
        if m := _S.match(line):
            values = parse_pairs(m.group(1))
            states[int(values["iteration"])] = values
        if m := _Q.match(line):
            if final is not None:
                raise ValueError("duplicate execution summary")
            final = (m.group(1), int(m.group(2)))
    if initial is None or not commits or final is None:
        raise ValueError("missing frozen-system, physical commit or execution summary")
    if sorted(commits) != list(range(1, max(commits) + 1)):
        raise ValueError("committed iteration gap")
    if max(commits) != final[1]:
        raise ValueError("uncommitted final iteration")
    if final[0] != "NO":
        raise ValueError("candidate is not a validated converged BFS experiment")
    if any(i not in commits for i in states):
        raise ValueError("field period sample without physical commit")
    if triggers and any(i not in commits for i in triggers):
        raise ValueError("phase predictor used without committing a full trial")
    if triggers and min(triggers) >= 2:
        previous=states.get(min(triggers)-1)
        if not previous or previous.get("period_four_candidate") != "1":
            raise ValueError("phase predictor lacks a preceding committed orbit")
    return {"initial": initial, "commits": commits, "trigger_iterations": triggers,
            "final": final, "field_samples": len(states)}


def compare(baseline_text: str, trial_text: str, require_trigger: bool = False) -> dict:
    baseline=extract(baseline_text)
    trial=extract(trial_text)
    if baseline["initial"] != trial["initial"]:
        raise ValueError("initial coupled A,b,x0 fingerprint differs")
    if baseline["trigger_iterations"]:
        raise ValueError("baseline used full V44 predictor")
    triggers=trial["trigger_iterations"]
    if require_trigger and not triggers:
        raise ValueError("V44 full-field predictor never triggered")
    first=min(triggers) if triggers else None
    common=sorted(set(baseline["commits"]) & set(trial["commits"]))
    if not common or common != list(range(1, max(common)+1)):
        raise ValueError("unmatched committed-state history")
    # No V44 intervention should change physical results at earlier states.
    if first:
        for i in common:
            if i>=first: break
            if (abs(trial["commits"][i]-baseline["commits"][i]) >
                    1e-10*max(1.,abs(baseline["commits"][i]))):
                raise ValueError("physical trajectory drifted BEFORE first V44 trial")
    last=common[-1]
    data={"status":"UNQUALIFIED_EXPERIMENTAL_PREDICTOR",
          "physical_converged":False,
          "initial_A_b_x0":baseline["initial"],
          "common_committed_states":len(common),
          "first_full_sst_predictor":first,
          "v44_predictor_attempts":len(triggers),
          "common_last_iteration":last,
          "baseline_common_last_merit":baseline["commits"][last],
          "trial_common_last_merit":trial["commits"][last],
          "baseline_final_iteration":baseline["final"][1],
          "trial_final_iteration":trial["final"][1]}
    if first and first in baseline["commits"]:
        data["baseline_first_trial_merit"]=baseline["commits"][first]
        data["v44_first_trial_merit"]=trial["commits"][first]
        data["merit_ratio_at_first_trial"]=(trial["commits"][first] /
                                            max(baseline["commits"][first],1e-300))
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline",type=Path,required=True)
    parser.add_argument("--trial",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--require-trigger",action="store_true")
    args=parser.parse_args()
    result=compare(args.baseline.read_text(),args.trial.read_text(),args.require_trigger)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(f"V44_COMPARE {result['status']} common={result['common_committed_states']} "
          f"trigger={result['first_full_sst_predictor']}")

if __name__ == "__main__":
    main()
