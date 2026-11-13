"""Fail-closed comparison of the SAME frozen 6N X and the real Jv(X).

Five Schur variants have a native 4N inverse available. The sixth (PCD)
failed its 4N linear solve and is explicitly marked unavailable rather than
quietly substituting another preconditioner. This is NOT a convergence gate.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
import sys

NAMES=("block_local", "lsc", "bfbt", "simple", "simplec")
RECORD=re.compile(r"^FPHYS_JFNK_V52 [^\n]+$",re.M)
SOURCE=re.compile(r"^COUPLED_MATRIX_FINGERPRINT_V40 [^\n]+$",re.M)

def fields(line: str) -> dict[str,str]:
    return dict(pair.split("=",1) for pair in line.split() if "=" in pair)

def audit_logs(logs: dict[str,str], audit_4n: dict) -> dict:
    if (audit_4n.get("status")!="COMPARABLE_4N_ONLY" or
        audit_4n.get("all_six_same_original_4N_A_b_x0_SHA256") is not True or
        audit_4n.get("n_cells")!=152):
        raise ValueError("missing verified common 4N A,b,x0 capture")
    if set(logs)!=set(NAMES):
        raise ValueError("five native 4N Schur logs required, PCD 4N is unavailable")
    entries=[]
    for name in NAMES:
        log=logs[name]
        found=RECORD.findall(log)
        sources=SOURCE.findall(log)
        if len(found)!=1 or len(sources)!=1:
            raise ValueError(f"incomplete or duplicate physical JFNK evidence for {name}")
        parts=fields(found[0]); base=fields(sources[0])
        if (parts.get("selected_schur")!=name or
            parts.get("frozen_prestep")!="1" or
            parts.get("native_schur_ready")!="1" or
            parts.get("production_committed")!="0" or
            parts.get("qualified")!="0" or
            parts.get("step")!="0" or
            not parts.get("common_state_fingerprint","").isdigit()):
            raise ValueError(f"incorrect provenance or false qualification: {name}")
        try:
            relative=float(parts["linear_true_relative"])
            merit=float(parts["merit_before"])
            krylov_evaluations=int(parts["evaluations"])
        except (ValueError,KeyError) as exc:
            raise ValueError(f"invalid JFNK certificate fields: {name}") from exc
        if not (math.isfinite(relative) and relative>=0 and math.isfinite(merit)
                and merit>=0 and krylov_evaluations>=2):
            raise ValueError(f"invalid nonfinite JFNK result: {name}")
        entries.append({"schur":name,"initial_6N_state_fingerprint":parts["common_state_fingerprint"],
                        "jv_true_linear_relative":relative,"physical_merit_before":merit,
                        "private_trial_accepted":parts.get("private_accepted")=="1",
                        "frozen_matrix_fingerprint":tuple(base.get(f) for f in ("A","b","x0")),
                        "candidate_evaluations":krylov_evaluations})
    if len({e["initial_6N_state_fingerprint"] for e in entries})!=1:
        raise ValueError("six-field initial state differs across Schur methods")
    if len({e["frozen_matrix_fingerprint"] for e in entries})!=1:
        raise ValueError("coupled operator fingerprints differ across Schur methods")
    if max(e["physical_merit_before"] for e in entries)-min(e["physical_merit_before"] for e in entries)>1e-10:
        raise ValueError("physical residual merit disagrees for identical frozen X")
    if any(e["private_trial_accepted"] and e["jv_true_linear_relative"]>0.3 for e in entries):
        raise ValueError("unqualified Krylov direction accepted")
    for entry in entries:entry.pop("frozen_matrix_fingerprint")
    return {"schema":"CFDX-V52-FROZEN-COMMON-6N-JFNK/1",
            "same_initial_6N_state":True,"common_state_fingerprint":entries[0]["initial_6N_state_fingerprint"],
            "same_original_4N_A_b_x0":True,
            "same_matrix_free_6N_F_definition":True,
            "explicit_6N_Jacobian_matrix_captured":False,
            "PCD_6N_not_testable_due_to_4N_failure":True,
            "all_six_6N_methods_qualified":False,
            "physically_qualified":False,
            "status":"SAME_FROZEN_X_5_METHODS_NO_NEWTON_CONVERGENCE"
                if all(not e["private_trial_accepted"] for e in entries)
                else "SAME_FROZEN_X_EXPERIMENTAL_PRIVATE_STEPS",
            "methods":[{"schur":"pcd","6n_jfnk":"UNAVAILABLE_PCD_4N_FAIL"}]+entries}

def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory",type=Path,required=True)
    parser.add_argument("--audit-4n",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    opts=parser.parse_args()
    try:
        logs={name:(opts.directory/f"frozen_{name}.log").read_text(encoding="utf-8") for name in NAMES}
        previous=json.loads(opts.audit_4n.read_text(encoding="utf-8"))
        report=audit_logs(logs,previous)
        opts.output.parent.mkdir(parents=True,exist_ok=True)
        opts.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
        print(report["status"])
        return 0
    except (ValueError,OSError,KeyError,json.JSONDecodeError) as err:
        print(f"V52 frozen six-field audit failure: {err}",file=sys.stderr)
        return 2

if __name__=="__main__":raise SystemExit(main())
