"""Strict non-regression / falsification audit for private V71 BFS152 trials.

Read-only experiments can never certify a production convergence. The audit
compares V65/V69 production fingerprints line by line with a supplied baseline.
"""
import argparse
import json
import math
from pathlib import Path

HISTORICAL = (
    "GLOBAL_SAME_OPERATOR_TRIAL_V65",
    "GLOBAL_SAME_OPERATOR_COMMITTED_V65",
    "BFS_SST_RC_PICARD_V69",
    "INCOMPRESSIBLE_ITER",
)
PRIVATE = ("BFS_JOINT_PICARD_V71", "BFS_JOINT_FACE3_V71")


def extract(log, marker):
    return [line.strip() for line in log.splitlines() if line.startswith(marker + " ")]


def fields(line):
    out = {}
    for chunk in line.split()[1:]:
        if "=" in chunk:
            key, value = chunk.split("=", 1)
            out[key] = value
    return out


def audit(baseline, candidate):
    comparison = {}
    for marker in HISTORICAL:
        old = extract(baseline, marker)
        new = extract(candidate, marker)
        if old != new:
            raise ValueError(f"V71 changed historical marker {marker}: {len(old)} vs {len(new)}")
        comparison[marker] = len(new)
    attempts = extract(candidate, PRIVATE[0])
    faces = extract(candidate, PRIVATE[1])
    if not attempts:
        raise ValueError("V71 joint solver was not executed")
    allowed_cells = {"0", "12"}
    for line in attempts:
        data = fields(line)
        if data.get("committed") != "0" or data.get("qualified") != "0" or data.get("private_trial") != "1":
            raise ValueError("V71 private audit illegally published or qualified a candidate")
        if data.get("iteration") != "23" or data.get("retry") != "8":
            raise ValueError("V71 evaluated an unapproved BFS152 iteration")
        if "rejected" not in data:
            for key in ("mass_before", "mass_after_predictor", "mass_after_frozen_pressure",
                        "mass_after_sst", "merit_ratio", "worst_growth", "omega_after"):
                if key not in data or not math.isfinite(float(data[key])):
                    raise ValueError("V71 incomplete/nonfinite physical merit or mass audit: " + key)
            if data.get("physically_eligible") == "1" and (
                data.get("physical_mass_closed") != "1" or data.get("original_guards_ok") != "1"
            ):
                raise ValueError("V71 falsely marked a candidate eligible")
    face_groups = {}
    maximum_budget_closure = 0.0
    for line in faces:
        data = fields(line)
        if data.get("cell") not in allowed_cells:
            raise ValueError("V71 unexpected face3 SST cell")
        if data.get("committed") != "0" or data.get("qualified") != "0":
            raise ValueError("V71 face audit falsely published a solution")
        for key in ("phi3_before", "phi3_after", "rAUx_before", "rAUx_after",
                    "nut_before", "nut_after", "omega_raw_before", "omega_raw_after",
                    "convection_delta", "diffusion_delta", "production_delta",
                    "destruction_delta", "cross_delta"):
            if key not in data or not math.isfinite(float(data[key])):
                raise ValueError("V71 missing/nonfinite face3 source term: " + key)
        if float(data["rAUx_before"]) <= 0 or float(data["rAUx_after"]) <= 0 or (float(data["nut_before"]) < 0 or float(data["nut_after"]) < 0):
            raise ValueError("V71 invalid face3 rAU/viscosity")
        omega_delta=float(data["omega_raw_after"])-float(data["omega_raw_before"])
        budget_delta=sum(float(data[key]) for key in
                         ("convection_delta", "diffusion_delta", "production_delta",
                          "destruction_delta", "cross_delta"))
        defect=abs(omega_delta-budget_delta)
        maximum_budget_closure=max(maximum_budget_closure,defect)
        if defect > 1e-8*max(1.,abs(omega_delta),abs(budget_delta)):
            raise ValueError("V71 independent omega source budget did not close")
        key=(data.get("gamma"),data.get("pressure_fraction"))
        face_groups.setdefault(key,[]).append(data["cell"])
    for cells in face_groups.values():
        if sorted(cells) != ["0", "12"]:
            raise ValueError("V71 incomplete cell0/cell12 face3 pair")
    if len(faces) != 2*sum("rejected" not in fields(line) for line in attempts):
        raise ValueError("V71 face3 rows do not match completed probes")
    successful=[fields(line) for line in attempts if "rejected" not in fields(line)]
    best=min(successful,key=lambda row:float(row["merit_ratio"])) if successful else None
    guarded=[row for row in successful if row.get("original_guards_ok")=="1"]
    best_guarded=min(guarded,key=lambda row:float(row["mass_after_sst"])) if guarded else None
    return {"status": "diagnostic_only_no_convergence_certificate", "historical_identity": comparison,
            "joint_trial_attempts": len(attempts), "successful_joint_trial_evaluations":
            sum("rejected" not in fields(line) for line in attempts),
            "face_cell_budget_rows": len(faces),
            "maximum_signed_omega_budget_closure": maximum_budget_closure,
            "best_merit_trial": {k: best.get(k) for k in
                ("gamma", "pressure_fraction", "merit_ratio", "worst_growth",
                 "mass_after_sst", "omega_after", "physically_eligible")} if best else None,
            "best_guarded_trial": {k: best_guarded.get(k) for k in
                ("gamma", "pressure_fraction", "merit_ratio", "worst_growth",
                 "mass_after_sst", "omega_after", "physically_eligible")} if best_guarded else None,
            "candidates_passing_v65_guards": sum(row.get("original_guards_ok")=="1" for row in successful),
            "candidates_conservative_and_guarded": sum(row.get("physically_eligible")=="1" for row in successful),
            "accepted_new_production_states": 0,
            "physical_convergence_qualified": False}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--log", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    result = audit(args.baseline.read_text(errors="replace"),
                   args.log.read_text(errors="replace"))
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()