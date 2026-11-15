#!/usr/bin/env python3
"""Fail-closed parser for the four-corner BFS V64 read-only physical diagnostic.

Requires all six signed RAW residual blocks and independently closed RC at all
four corners. A certificate from this script is DIAGNOSTIC ONLY, never a
nonlinear convergence gate or a claim that isolated frozen solvers converged.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path

SUMMARY = "BFS_CAUSAL_V64"
BLOCK = "BFS_CAUSAL_BLOCK_V64"
REQUIRED_SUMMARY = {"iteration", "retry", "cells", "flow_changed", "sst_changed",
    "relative_flow_delta", "relative_sst_delta", "recomposition_error",
    "rc00", "rc10", "rc01", "rc11", "qualified", "committed",
    "four_corners_reassembled", "phi_flow", "phi_sst",
    "phi_interaction", "phi_interaction_face", "nut_flow", "nut_sst",
    "nut_interaction", "nut_interaction_cell"}
REQUIRED_BLOCK = {"iteration", "block", "raw00", "raw10", "raw01",
    "raw11", "flow_only_delta", "sst_only_delta", "interaction",
    "interaction_relative", "worst_after_cell", "worst_interaction_cell",
    "qualified", "committed"}


def tokens(line: str) -> dict[str, str]:
    return dict(re.findall(r"([A-Za-z][A-Za-z_0-9]*)=([^\s]+)", line))


def parse(log: str) -> dict:
    summaries, blocks = [], {}
    for line in log.splitlines():
        if line.startswith(SUMMARY + " "):
            kv = tokens(line)
            if REQUIRED_SUMMARY - kv.keys():
                raise ValueError("incomplete V64 causal summary")
            summaries.append(kv)
        elif line.startswith(BLOCK + " "):
            kv = tokens(line)
            if REQUIRED_BLOCK - kv.keys():
                raise ValueError("incomplete V64 causal block")
            index = (int(kv["iteration"]), int(kv["block"]))
            if index in blocks:
                raise ValueError("duplicate V64 causal residual block")
            blocks[index] = kv
    if not summaries:
        raise ValueError("no V64 causal records: not executed")
    if len({int(s["iteration"]) for s in summaries}) != len(summaries):
        raise ValueError("duplicate V64 causal summaries")
    reports = []
    for s in summaries:
        it = int(s["iteration"])
        n = int(s["cells"])
        if n < 1 or n > 152 or int(s["retry"]) != 0:
            raise ValueError("V64 iteration/mesh/prestate is outside bounded protocol")
        if (s["qualified"] != "0" or s["committed"] != "0" or
            s["four_corners_reassembled"] != "1"):
            raise ValueError("V64 records falsely claim committed or certified state")
        def number(kv, key):
            x = float(kv[key])
            if not math.isfinite(x) or x < 0.0:
                raise ValueError("non-finite or negative V64 metric: " + key)
            return x
        for key in ("phi_flow", "phi_sst", "phi_interaction",
                    "nut_flow", "nut_sst", "nut_interaction"):
            number(s,key)
        for key in ("phi_interaction_face", "nut_interaction_cell"):
            if int(s[key])<0 or (key=="nut_interaction_cell" and int(s[key])>=n):
                raise ValueError("V64 operator hotspot outside valid range")
        for key in ("rc00", "rc10", "rc01", "rc11"):
            if number(s,key) > 1e-10:
                raise ValueError("V64 RC fixed point not closed at one corner")
        for key in ("flow_changed", "sst_changed"):
            if s[key] not in {"0", "1"}:
                raise ValueError("V64 change flag invalid")
        flow_delta = number(s,"relative_flow_delta")
        sst_delta = number(s,"relative_sst_delta")
        if (s["flow_changed"] == "0") != (flow_delta == 0.0):
            raise ValueError("V64 flow change contradicts its magnitude")
        if (s["sst_changed"] == "0") != (sst_delta == 0.0):
            raise ValueError("V64 SST change contradicts its magnitude")
        rows = []
        for b in range(6):
            kv=blocks.get((it,b))
            if kv is None:
                raise ValueError("missing V64 physical block " + str(b))
            if kv["qualified"] != "0" or kv["committed"] != "0":
                raise ValueError("V64 block falsely claims acceptance")
            for key in ("raw00","raw10","raw01","raw11","flow_only_delta",
                        "sst_only_delta","interaction","interaction_relative"):
                number(kv,key)
            for key in ("worst_after_cell","worst_interaction_cell"):
                if not 0 <= int(kv[key]) < n:
                    raise ValueError("V64 hotspot cell outside mesh")
            rows.append({"block":b,**{key:float(kv[key]) for key in
                ("raw00","raw10","raw01","raw11","flow_only_delta",
                 "sst_only_delta","interaction","interaction_relative")},
                "hotspot_after":int(kv["worst_after_cell"]),
                "hotspot_interaction":int(kv["worst_interaction_cell"])})
        if number(s,"recomposition_error") > 1e-7 * max(
            1.0,*(max(row[k] for k in ("raw00","raw10","raw01","raw11"))
                  for row in rows)):
            raise ValueError("V64 four-corner closure not reproducible")
        reports.append({"iteration":it,"cells":n,"flow_changed":s["flow_changed"]=="1",
            "sst_changed":s["sst_changed"]=="1", "blocks":rows,
            "max_rc_defect":max(float(s[key]) for key in
                ("rc00","rc10","rc01","rc11")),
            "rhie_chow_coupling":{key:float(s[key]) for key in
                ("phi_flow","phi_sst","phi_interaction")},
            "eddy_viscosity_coupling":{key:float(s[key]) for key in
                ("nut_flow","nut_sst","nut_interaction")}})
    if len(blocks) != 6 * len(reports):
        raise ValueError("orphan V64 physical blocks")
    return {"schema":"CFDX-BFS-CAUSAL-COUNTERFACTUAL/V64",
        "status":"DIAGNOSTIC_ONLY_NOT_A_FROZEN_SOLVE",
        "newton_accepted":False,"physical_qualification":False,
        "iterations":reports}


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=parse(args.log.read_text(encoding="utf-8",errors="replace"))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("V64_CAUSAL_EVIDENCE_COMPLETE_DIAGNOSTIC_ONLY")
    return 0

if __name__=="__main__":
    raise SystemExit(main())