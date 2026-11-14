"""Fail-closed, first-attempt-only V24/V25 SST globalization diagnostic.

This parser does not certify a CFD solution. It separates a linear GMRES
candidate, dynamic-physical line search and the actual committed decision.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path


def first_attempt(text: str) -> dict:
    data = text.split("CFDX nonlinear iteration 1 rejected", 1)[0]
    jfnk = [line for line in data.splitlines()
            if line.startswith("SST_JFNK_V22 iteration=10 ")]
    physical = [line for line in data.splitlines()
                if line.startswith("SST_TRANSPORT_DIAGNOSTIC iteration=10 ")]
    gates = [line for line in data.splitlines()
             if line.startswith("SST_ACCEPTANCE_GATES iteration=10 ")]
    if len(jfnk) != 1 or len(physical) != 1 or len(gates) != 1:
        raise ValueError("V25 diagnosis requires one complete first-attempt SST iteration 10")

    def parse(line: str) -> dict:
        return dict(re.findall(r"([A-Za-z][A-Za-z_0-9]*)=([^\s]+)", line))

    j, p, g = map(parse, (jfnk[0], physical[0], gates[0]))
    need = ("spatial_linear_true_relative", "spatial_probe", "selected",
            "old_merit", "new_merit", "v25_active", "v25_direction",
            "v25_accepted", "v25_alpha", "v25_best_trial_merit",
            "v25_wall_corrected")
    if any(k not in j for k in need):
        raise ValueError("incomplete globalizer evidence")
    num_keys = ("spatial_linear_true_relative", "old_merit", "new_merit",
                "v25_alpha", "v25_best_trial_merit")
    for key in num_keys:
        j[key] = float(j[key])
    for key in ("v25_active", "v25_direction", "v25_accepted",
                "v25_wall_corrected", "spatial_probe"):
        j[key] = int(j[key])
    for val in (j[k] for k in num_keys):
        if not math.isfinite(val):
            raise ValueError("nonfinite measured merit or linear residual")
    if j["v25_active"] != 1 or j["spatial_probe"] != 1:
        raise ValueError("V25 was not actually enabled on the first attempt")
    if j["v25_accepted"] and not (j["v25_direction"] and
                                   j["v25_alpha"] > 0 and
                                   j["v25_best_trial_merit"] < j["old_merit"]):
        raise ValueError("V25 claimed a non-descending log step")
    if j["selected"] == "jfnk_spatial_log_v25" and not j["v25_accepted"]:
        raise ValueError("unavailable V25 log step selected")
    if j["new_merit"] > j["old_merit"]:
        raise ValueError("physical step increased merit")
    if "tolerance" not in g or "qualified_bfs" not in g:
        raise ValueError("missing unchanged physical convergence gate")
    if float(g["tolerance"]) != 1e-8 or int(g["qualified_bfs"]) != 0:
        raise ValueError("unexpected tolerance or false physical qualification")
    if "k_local_relative" not in p or "omega_local_relative" not in p:
        raise ValueError("missing reassembled physical residual")
    for key in ("k_local_relative", "omega_local_relative"):
        if not math.isfinite(float(p[key])):
            raise ValueError("nonfinite physical residual")
    return {
        "iteration": 10,
        "spatial_krylov_relative": j["spatial_linear_true_relative"],
        "log_direction_available": bool(j["v25_direction"]),
        "v25_physical_log_step_found": bool(j["v25_accepted"]),
        "v25_alpha": j["v25_alpha"],
        "v25_log_merit": j["v25_best_trial_merit"],
        "old_merit": j["old_merit"],
        "selected_merit": j["new_merit"],
        "selected_direction": j["selected"],
        "wall_blocks_corrected": j["v25_wall_corrected"],
        "k_local_relative": float(p["k_local_relative"]),
        "omega_local_relative": float(p["omega_local_relative"]),
        "sst_full_convergence_certified": False,
        "bfs_global_convergence_certified": False,
        "openfoam_equivalence_certified": False,
        "nasa_physical_validation_certified": False,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log",type=Path)
    parser.add_argument("--output",type=Path)
    a=parser.parse_args()
    result=first_attempt(a.log.read_text(errors="replace"))
    rendered=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
    if a.output:
        a.output.write_text(rendered)
    print(rendered,end="")


if __name__=="__main__": main()
