"""V74 fail-closed verifier for same-state fixed-direction and wall evidence.

Running with no arguments executes adversarial synthetic regressions.
With --baseline/--log it compares authoritative V65 trial/certificates and
rejects any leaked commit or inconsistent mass/merit record.
"""
import argparse
import json
import math
import re
from pathlib import Path


def parse_kv(line):
    return dict(re.findall(r"([a-zA-Z][a-zA-Z0-9_]*)=([^\s]+)", line))


def lines_with(text, token):
    return [line for line in text.splitlines() if line.startswith(token)]


def audit(baseline, log, require_wall=True):
    baseline_tokens=("GLOBAL_SAME_OPERATOR_TRIAL_V65 ",
                     "GLOBAL_SAME_OPERATOR_COMMITTED_V65 ",
                     "INCOMPRESSIBLE_ITER ")
    for tag in baseline_tokens:
        old=lines_with(baseline,tag)
        new=lines_with(log,tag)
        if not old or old!=new:
            raise ValueError("V74 changes authoritative production history: "+tag)
    attempts=lines_with(log,"BFS_FIXED_DIRECTION_V74 ")
    summary=lines_with(log,"BFS_FIXED_DIRECTION_SUMMARY_V74 ")
    if not attempts or len(summary)!=1:
        raise ValueError("V74 missing bounded continuation evidence")
    if len(attempts)>8:
        raise ValueError("V74 exceeded single direction budget")
    selected=0
    previous_mass=None
    for index,line in enumerate(attempts,1):
        a=parse_kv(line)
        for key in ("iteration","cycle","evaluated","selected","mass_before",
                    "mass_after","merit_before","merit_after","guards",
                    "mass_descent","relevant_groups_nonincreasing","phi3_before",
                    "phi3_after","rAU0_before","rAU0_after","rAU12_before",
                    "rAU12_after","nut0_before","nut0_after","nut12_before",
                    "nut12_after","mass_closed","six_closed","committed","qualified"):
            if key not in a:
                raise ValueError("V74 incomplete row "+key)
        if int(a["iteration"])!=23 or int(a["cycle"])!=index or \
           a["committed"]!="0" or a["qualified"]!="0":
            raise ValueError("V74 false iteration, cycle, or qualification")
        vals={key:float(a[key]) for key in ("mass_before","mass_after",
            "merit_before","merit_after","phi3_before","phi3_after",
            "rAU0_before","rAU0_after","rAU12_before","rAU12_after",
            "nut0_before","nut0_after","nut12_before","nut12_after")}
        if not all(math.isfinite(v) for v in vals.values()):
            raise ValueError("V74 nonfinite diagnostics")
        if previous_mass is not None and abs(vals["mass_before"]-previous_mass)>1e-10:
            raise ValueError("V74 stale intermediate input")
        if a["selected"]=="1":
            selected+=1
            if not (a["evaluated"]==a["guards"]==a["mass_descent"]==
                    a["relevant_groups_nonincreasing"]=="1"):
                raise ValueError("V74 private state selected without every guard")
            if not (vals["mass_after"]<vals["mass_before"] or
                    vals["mass_after"]<=1e-10):
                raise ValueError("V74 selected nonconservative step")
            if not vals["merit_after"]<vals["merit_before"]:
                raise ValueError("V74 selected a non-descending merit")
            previous_mass=vals["mass_after"]
        else:
            if index!=len(attempts):
                raise ValueError("V74 continued after a rejected private state")
            previous_mass=vals["mass_before"]
    s=parse_kv(summary[0])
    if s.get("committed")!="0" or s.get("qualified")!="0" or \
       s.get("native_momentum")!="0" or s.get("native_sst")!="0" or \
       int(s.get("private_accepted","-1"))!=selected or \
       int(s.get("cycles","-1"))!=len(attempts):
        raise ValueError("V74 fake native certificate or summary mismatch")
    if selected>0 and abs(float(s["mass_after"])-previous_mass)>1e-10:
        raise ValueError("V74 final mass certificate differs from accepted state")
    walls=lines_with(log,"BFS_WALL_REGIME_V74 ")
    if require_wall:
        if len(walls)!=1:
            raise ValueError("V74 missing or duplicated wall regime certificate")
        w=parse_kv(walls[0]);count=int(w["faces"])
        if count<1 or sum(int(w[key]) for key in
              ("wall_resolved","buffer","log_layer","out_of_range"))!=count:
            raise ValueError("V74 missing wall faces from yplus classification")
        if w.get("openfoam_equivalent")!="0" or w.get("qualified")!="0":
            raise ValueError("V74 claims false OpenFOAM equivalence")
        for key in ("min_yplus_spalding","max_yplus_spalding","omega_step_gap",
                    "omega_blended_gap","wall_nut_gap"):
            if not math.isfinite(float(w[key])) or float(w[key])<0:
                raise ValueError("V74 invalid yplus / omega / nut wall evidence")
    return {"authoritative_trial_lines":len(lines_with(log,baseline_tokens[0])),
            "private_attempts":len(attempts),"private_selected":selected,
            "wall_faces":int(parse_kv(walls[0])["faces"]) if walls else 0,
            "production_commits":0,"convergence_certified":False}


def selftest():
    baseline=("GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=23 accepted=0\n"
              "GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=22 qualified=0\n"
              "INCOMPRESSIBLE_ITER iter=23 momentum=0.1\n")
    row=("BFS_FIXED_DIRECTION_V74 iteration=23 cycle=1 evaluated=1 selected=1 "
         "mass_before=0.0006 mass_after=0.0005 merit_before=0.2 merit_after=0.19 "
         "guards=1 mass_descent=1 relevant_groups_nonincreasing=1 "
         "phi3_before=0.002 phi3_after=0.001 rAU0_before=1 rAU0_after=1 "
         "rAU12_before=1 rAU12_after=1 nut0_before=0.1 nut0_after=0.1 "
         "nut12_before=0.1 nut12_after=0.1 mass_closed=0 six_closed=0 "
         "committed=0 qualified=0 reason=none\n")
    summary=("BFS_FIXED_DIRECTION_SUMMARY_V74 iteration=23 cycles=1 "
             "private_accepted=1 mass_before=0.0006 mass_after=0.0005 "
             "mass_closed=0 native_momentum=0 native_sst=0 "
             "committed=0 qualified=0 status=bounded_cycle_budget_exhausted\n")
    wall=("BFS_WALL_REGIME_V74 iteration=23 faces=2 wall_resolved=0 buffer=1 "
          "log_layer=1 out_of_range=0 min_yplus_spalding=10 "
          "max_yplus_spalding=100 omega_step_gap=0.3 omega_blended_gap=0.2 "
          "wall_nut_gap=0.4 openfoam_equivalent=0 qualified=0\n")
    valid=baseline+row+summary+wall
    assert audit(baseline,valid)["private_selected"]==1
    alterations=[valid.replace("selected=1","selected=0"),
                 valid.replace("mass_after=0.0005","mass_after=0.0007"),
                 valid.replace("guards=1","guards=0"),
                 valid.replace("relevant_groups_nonincreasing=1","relevant_groups_nonincreasing=0"),
                 valid.replace("native_momentum=0","native_momentum=1"),
                 valid.replace("openfoam_equivalent=0","openfoam_equivalent=1"),
                 valid.replace("wall_resolved=0","wall_resolved=1"),
                 valid.replace("GLOBAL_SAME_OPERATOR_TRIAL_V65","CHANGED_TRIAL"),
                 valid.replace("omega_step_gap=0.3","omega_step_gap=nan"),
                 valid.replace("committed=0","committed=1"),
                 valid.replace("cycle=1","cycle=2")]
    for bad in alterations:
        try:audit(baseline,bad)
        except (ValueError,KeyError):continue
        raise AssertionError("V74 falsification was accepted")
    print("V74_EVIDENCE_FALSIFICATION_PASS cases=11")


if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument('--baseline');a.add_argument('--log');a.add_argument('--output')
    opts=a.parse_args()
    if opts.baseline or opts.log:
        if not opts.baseline or not opts.log:
            a.error('--baseline and --log are both required')
        result=audit(Path(opts.baseline).read_text(),Path(opts.log).read_text())
        if opts.output:Path(opts.output).write_text(json.dumps(result,indent=2))
        print(json.dumps(result,sort_keys=True))
    else:selftest()
