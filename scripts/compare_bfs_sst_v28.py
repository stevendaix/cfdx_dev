#!/usr/bin/env python3
"""Compare only the first attempted BFS outer state. Never a V&V certificate."""
from __future__ import annotations
import argparse,json,math,sys
from pathlib import Path
from diagnose_bfs_sst_v26 import kv

_ALLOWED={"picard","jfnk","jfnk_block","jfnk_spatial",
          "jfnk_spatial_log_v25","jfnk_smooth_v26",
          "sst_subspace_trust_v27","sst_adaptive_trust_v28","none"}

def trajectory(contents: str, expect_v28: bool) -> dict:
    marker="CFDX nonlinear iteration 1 rejected"
    if marker not in contents:
        raise ValueError("first attempted outer state has no explicit rejection marker")
    first=contents.split(marker,1)[0]
    raw=[kv(line) for line in first.splitlines()
         if "SST_JFNK_V22 iteration=" in line]
    if not raw:
        raise ValueError("missing per-iteration JFNK physics trace")
    trace=[]
    for row in raw:
        try:
            number=int(row["iteration"])
            old=float(row["old_merit"])
            new=float(row["new_merit"])
            active=(row["v28_active"]=="1")
            selected=row["selected"]
            committed=(row["v28_committed"]=="1")
        except (ValueError,KeyError) as e:
            raise ValueError("invalid mandatory numerical trace") from e
        if active!=expect_v28 or selected not in _ALLOWED or not number>0:
            raise ValueError("unexpected experiment branch or iteration")
        if not (math.isfinite(old) and math.isfinite(new) and
                old>=0 and new>=0 and
                (selected=="none" and abs(new-old)<=1e-8*max(1.,old) or
                 selected!="none" and new<old)):
            raise ValueError("invalid final physical FV acceptance")
        if committed and (selected not in {"sst_subspace_trust_v27",
                                          "sst_adaptive_trust_v28"} or
                          row.get("v28_commit_source") not in {"v27","v28"}):
            raise ValueError("radius updated without a committed trust candidate")
        if trace and number<=trace[-1]["iteration"]:
            raise ValueError("nonmonotonic or duplicated first-attempt trace")
        trace.append({"iteration":number,"before":old,"after":new,
                      "selected":selected,"v28_rescue_probe":
                            row.get("v28_rescue_probe")=="1",
                      "v28_committed":committed})
    return {"first_attempt_rejected":True,"outer_accepted":False,
            "sst_converged":False,"bfs_qualified":False,
            "v28_active":expect_v28,"last_logged_linf":trace[-1]["after"],
            "last_logged_iteration":trace[-1]["iteration"],"trace":trace}

def compare(off: str,on: str) -> dict:
    a,b=trajectory(off,False),trajectory(on,True)
    a_by={row["iteration"]:row for row in a["trace"]}
    b_by={row["iteration"]:row for row in b["trace"]}
    for i in (1,10,20):
        if i not in a_by or i not in b_by:
            raise ValueError(f"missing baseline check iteration {i}")
        aa,bb=a_by[i],b_by[i]
        if aa["selected"]!=bb["selected"] or not math.isclose(
            aa["after"],bb["after"],rel_tol=1e-11,abs_tol=1e-12):
            raise ValueError(f"pre-rescue physical mismatch at iteration {i}")
    return {"baseline":a,"adaptive":b,
            "extra_linf_contraction_after_stall":
                a["last_logged_linf"]-b["last_logged_linf"],
            "converged":False,"qualified":False,
            "note":"SST inner iteration remained non-converged; outer solve rolled back"}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("off",type=Path);p.add_argument("on",type=Path)
    p.add_argument("-o","--output",type=Path)
    args=p.parse_args()
    result=compare(args.off.read_text(),args.on.read_text())
    dump=json.dumps(result,indent=2,ensure_ascii=False)
    if args.output:args.output.write_text(dump+"\n")
    print(dump)
