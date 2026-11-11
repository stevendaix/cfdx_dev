#!/usr/bin/env python3
"""Fail-closed A/B proof for experimental SST branch preconditioning V34.

This tool never certifies BFS physics or a globally accepted outer iteration.
The reference HDF5 and V34 HDF5 must be byte-for-byte identical on geometry;
only one explicit numerics option may vary.
"""
import argparse
import json
import math
import re
from pathlib import Path

PAIR = re.compile(r"([a-zA-Z_][a-zA-Z_0-9]*)=([^\s]+)")
FLAG = "sst_branch_preconditioner_v34_experimental"


def first_outer_attempt(text):
    marker = "CFDX nonlinear iteration 1 rejected"
    if marker not in text:
        raise ValueError("first outer rejection not proven; cannot certify comparison")
    return text.split(marker, 1)[0]


def sample(text, prefix, iteration="10"):
    matches = []
    for line in first_outer_attempt(text).splitlines():
        if line.startswith(prefix + " "):
            fields = dict(PAIR.findall(line))
            if fields.get("iteration") == iteration:
                matches.append(fields)
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {prefix} iteration={iteration}")
    return matches[0]


def verify_cases(reference, experimental):
    import h5py
    with h5py.File(reference) as a, h5py.File(experimental) as b:
        if set(a.keys()) != set(b.keys()):
            raise ValueError("HDF5 dataset lists differ")
        for key in a.keys():
            if a[key][()].tobytes() != b[key][()].tobytes():
                raise ValueError(f"HDF5 geometry or field differs: {key}")
        attrs_a, attrs_b = dict(a.attrs), dict(b.attrs)
        for key in attrs_a.keys() | attrs_b.keys():
            if key == "case_setup_json":
                continue
            if key not in attrs_a or key not in attrs_b or attrs_a[key] != attrs_b[key]:
                raise ValueError(f"HDF5 metadata differs: {key}")
        config_a = json.loads(attrs_a["case_setup_json"])
        config_b = json.loads(attrs_b["case_setup_json"])
    settings_a = config_a["numerics"]["raw_settings"]
    settings_b = config_b["numerics"]["raw_settings"]
    if settings_a.pop(FLAG, "false") != "false" or settings_b.pop(FLAG, None) != "true":
        raise ValueError("V34 flag not set as control/experiment")
    if config_a != config_b:
        raise ValueError("cases vary beyond V34 activation flag")


def _number(fields, key):
    number = float(fields[key])
    if not math.isfinite(number):
        raise ValueError(f"nonfinite {key}")
    return number


def compare(reference, experimental):
    ctrl = sample(reference, "SST_TRANSPORT_DIAGNOSTIC")
    trial = sample(experimental, "SST_TRANSPORT_DIAGNOSTIC")
    j0 = sample(reference, "SST_JFNK_V22")
    j1 = sample(experimental, "SST_JFNK_V22")
    g0 = sample(reference, "SST_ACCEPTANCE_GATES")
    g1 = sample(experimental, "SST_ACCEPTANCE_GATES")
    correction = sample(experimental, "SST_BRANCH_PRECOND_V34")
    for gate in (g0,g1):
        if gate.get("qualified_bfs") != "0" or gate.get("outer_update_ok") != "0":
            raise ValueError("unexpected qualified/accepted outer state")
        if _number(gate, "tolerance") != 1e-8:
            raise ValueError("physical SST tolerance changed")
    if correction.get("qualified") != "0" or int(correction["evaluated"]) > 50:
        raise ValueError("V34 correction falsely qualified or unbounded")
    if not correction.get("cells"):
        raise ValueError("V34 failed to sample a physical hotspot")
    selected0,selected1 = j0["selected"],j1["selected"]
    old_merit=max(_number(ctrl,"k_local_relative"),_number(ctrl,"omega_local_relative"))
    new_merit=max(_number(trial,"k_local_relative"),_number(trial,"omega_local_relative"))
    return {
        "status": "NOT_CONVERGED_FIRST_OUTER_REJECTED",
        "iteration":10,
        "same_physical_threshold":1e-8,
        "baseline_merit":old_merit,
        "experimental_merit":new_merit,
        "merit_difference":new_merit-old_merit,
        "selected_baseline":selected0,
        "selected_experimental":selected1,
        "cells":correction["cells"],
        "corrected_blocks":int(correction["corrected"]),
        "maximum_relative_difference":_number(correction,"max_difference"),
        "maximum_unilateral_gap":_number(correction,"max_unilateral_gap"),
        "source": "first_outer_attempt_only",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--control-log",required=True,type=Path)
    parser.add_argument("--experiment-log",required=True,type=Path)
    parser.add_argument("--control-hdf5",required=True,type=Path)
    parser.add_argument("--experiment-hdf5",required=True,type=Path)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    verify_cases(args.control_hdf5,args.experiment_hdf5)
    result=compare(args.control_log.read_text(),args.experiment_log.read_text())
    rendered=json.dumps(result,indent=2,sort_keys=True)
    if args.output:
        args.output.write_text(rendered+"\n")
    print(rendered)


if __name__ == "__main__":
    main()
