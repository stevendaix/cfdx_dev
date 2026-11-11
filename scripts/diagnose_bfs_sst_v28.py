"""Fail-closed diagnostic of the FIRST BFS outer attempt; never certifies BFS."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from diagnose_bfs_sst_v26 import kv

BRANCHES = {"picard", "jfnk", "jfnk_block", "jfnk_spatial",
            "jfnk_spatial_log_v25", "jfnk_smooth_v26",
            "sst_subspace_trust_v27", "sst_adaptive_trust_v28", "none"}

def first_attempt(log: str) -> dict:
    first = log.split("CFDX nonlinear iteration 1 rejected", 1)[0]
    traces = [kv(s) for s in first.splitlines() if "SST_JFNK_V22 iteration=10" in s]
    if not traces:
        raise ValueError("first attempt is missing its V28 iteration-10 trace")
    t = traces[0]
    if t.get("v28_active") != "1" or t.get("v27_active") != "1":
        raise ValueError("V28 opt-in absent")
    if t.get("qualified") != "0":
        raise ValueError("V28 diagnostics cannot certify the BFS")
    def num(key):
        try:
            v=float(t[key])
        except (ValueError,KeyError) as e:
            raise ValueError(f"missing or malformed V28 field {key}") from e
        if not math.isfinite(v):
            raise ValueError(f"nonfinite V28 field {key}")
        return v
    before,after = num("old_merit"),num("new_merit")
    selected=t.get("selected")
    if selected not in BRANCHES:
        raise ValueError("unknown selected SST branch")
    if selected != "none" and not after<before:
        raise ValueError("selected SST branch did not contract physical Linf")
    active=t.get("v28_accepted")=="1"
    committed=t.get("v28_committed")=="1"
    radius0,radius1=num("v28_radius_before"),num("v28_radius_after")
    if not (1/64<=radius0<=1.5 and 1/64<=radius1<=1.5):
        raise ValueError("V28 radius out of bounds")
    evaluations=num("v28_evaluations")
    candidates=num("v28_candidates")
    if not (evaluations>=0 and candidates>=0 and candidates<=evaluations):
        raise ValueError("bad V28 evaluation accounting")
    commit_source=t.get("v28_commit_source")
    if committed and (selected not in {"sst_adaptive_trust_v28",
                                  "sst_subspace_trust_v27"} or
                      num("alpha")!=1. or
                      commit_source not in {"v27","v28"} or
                      (commit_source=="v28" and
                       (not active or selected!="sst_adaptive_trust_v28")) or
                      (commit_source=="v27" and
                       (t.get("v27_accepted")!="1" or
                        selected!="sst_subspace_trust_v27" or
                        num("v27_ratio")<=1e-3))):
        raise ValueError("hypothetical V27/V28 proposal was counted as committed")
    if not committed and commit_source!="none":
        raise ValueError("false V28 trust commit source")
    if selected=="sst_adaptive_trust_v28" and not active:
        raise ValueError("unavailable V28 proposal selected")
    if active:
        l2=num("v28_l2")
        l2base=num("v26_l2_before")
        linf=num("v28_linf")
        ratio=num("v28_rho")
        weight=num("v28_jacobi_weight")
        step=num("v28_step")
        if not (0<step<=radius0*(1+1e-8) and 0<weight<1 and ratio>1e-3 and
                linf<before and l2<l2base):
            raise ValueError("V28 proposal does not meet its physical gates")
        if not committed and abs(radius0-radius1)>1e-12:
            raise ValueError("uncommitted V28 proposal changed trust radius")
    else:
        if selected=="sst_adaptive_trust_v28" or (committed and
                selected!="sst_subspace_trust_v27"):
            raise ValueError("invalid V28 selected or committed status")
    return {"v28_accepted":active,"v28_committed":committed,
            "selected":selected,"v28_evaluations":int(evaluations),
            "v28_candidates":int(candidates),"radius_before":radius0,
            "radius_after":radius1,"physical_linf_before":before,
            "physical_linf_after":after,"sst_converged":False,
            "bfs_converged":False,"nasa_validated":False}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('log',type=Path)
    a=p.parse_args()
    print(json.dumps(first_attempt(a.log.read_text()),indent=2))
