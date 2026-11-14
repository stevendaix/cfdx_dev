"""Fail-closed V53 6N Schur/SST audit. Algebraic reproducibility != BFS qualification.

Requires five actual V53 production runs plus the PCD failed attempt on
exactly the same 4N capture and the same frozen 6N initial state.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import pathlib
import re

METHODS = ("block_local", "lsc", "bfbt", "simple", "simplec", "pcd")

def parse(text: str, method: str) -> dict:
    def one(prefix):
        rows=[x for x in text.splitlines() if x.startswith(prefix)]
        if len(rows)!=1:
            raise ValueError(f"{method}: expected one {prefix}, got {len(rows)}")
        return dict(re.findall(r"([A-Za-z0-9_]+)=([^\s]+)", rows[0]))
    matrix=one("COUPLED_MATRIX_FINGERPRINT_V40 ")
    if matrix.get("rows")!="608" or matrix.get("nnz")!="7080":
        raise ValueError(f"{method}: wrong BFS152 captured 4N dimensions")
    if "Converged YES" in text:
        raise ValueError(f"{method}: falsely converged")
    if method!="pcd" and "Converged NO" not in text:
        raise ValueError(f"{method}: missing unconverged solver status")
    result={"method":method,"A":matrix.get("A"),"b":matrix.get("b"),
            "x0":matrix.get("x0"),"rows":608}
    for key in ("A","b","x0"):
        if not result[key]: raise ValueError(f"{method}: missing {key} fingerprint")
    if method=="pcd":
        if ("COUPLED_BLOCK_TRUE_RESIDUAL_V16" not in text or
            "coupled momentum-continuity solve did not converge" not in text):
            raise ValueError("PCD: 4N failure proof missing")
        if "FPHYS_JFNK_V52 " in text:
            raise ValueError("PCD: cannot promote failed 4N to a 6N setup")
        result.update({"available_for_6n":False,"linear_6n_certified":False})
        return result
    f=one("FPHYS_JFNK_V52 ")
    expected={"block_local":"block_local","lsc":"lsc","bfbt":"bfbt",
              "simple":"simple","simplec":"simplec"}
    for key,value in (("selected_schur",expected[method]),
                      ("v53_global_block","1"), ("native_schur_ready","1"),
                      ("frozen_prestep","1"), ("production_committed","0"),
                      ("qualified","0")):
        if f.get(key)!=value:
            raise ValueError(f"{method}: {key} does not match {value}")
    for field in ("linear_true_relative","arnoldi_relative","jv_linearity",
                  "merit_before","v53_locality_error"):
        if field not in f or not math.isfinite(float(f[field])) or float(f[field])<0:
            raise ValueError(f"{method}: invalid {field}")
    for key in ("v53_probes","v53_transport_nnz","common_state_fingerprint"):
        if not f.get(key,"").isdigit() or int(f[key])==0:
            raise ValueError(f"{method}: missing {key}")
    if int(f["v53_probes"])<=6 or int(f["v53_transport_nnz"])<=304:
        raise ValueError(f"{method}: missing real sparse SST block")
    if float(f["v53_locality_error"])>0.01:
        raise ValueError(f"{method}: colored local Jacobian contamination")
    if f.get("private_accepted")!="0":
        raise ValueError(f"{method}: unproved Newton acceptance")
    if float(f["linear_true_relative"])<=0.3:
        raise ValueError(f"{method}: unexpected success requires separate qualification")
    result.update({"available_for_6n":True,
        "linear_6n_certified":False,
        "common_state_fingerprint":f["common_state_fingerprint"],
        "linear_true_relative":float(f["linear_true_relative"]),
        "arnoldi_relative":float(f["arnoldi_relative"]),
        "jv_linearity_relative":float(f["jv_linearity"]),
        "transport_nnz":int(f["v53_transport_nnz"]),
        "p_mat_probes":int(f["v53_probes"])})
    return result

def audit(folder:pathlib.Path, captures:pathlib.Path|None=None) -> dict:
    records=[]
    for m in METHODS:
        p=folder/f"v53_{m}.log"
        if not p.is_file():raise ValueError(f"missing production log {p}")
        records.append(parse(p.read_text(errors="replace"),m))
    for col in ("A","b","x0"):
        if len({r[col] for r in records})!=1:
            raise ValueError(f"different {col} on six Schur")
    file_hashes={}
    if captures is not None:
        for tag in ("A","b","x0"):
            values={}
            for method in METHODS:
                file=captures/f"{method}_{tag}.mtx"
                if not file.is_file():
                    raise ValueError(f"missing 4N capture: {file}")
                values[method]=hashlib.sha256(file.read_bytes()).hexdigest()
            if len(set(values.values()))!=1:
                raise ValueError(f"different original 4N {tag} bytes")
            file_hashes[tag]=next(iter(values.values()))
    common=[r["common_state_fingerprint"] for r in records if r["available_for_6n"]]
    if len(set(common))!=1:
        raise ValueError("V53 Krylov baselines are different 6N states")
    return {"status":"INCOMPLETE_UNCONVERGED",
        "physically_qualified":False,
        "same_original_4n_A_b_x0":True,
        "same_initial_6n_candidate":True,
        "capture_sha256":file_hashes,
        "matrices_6n_explicitly_equal":False,
        "pcd_jfnk_result":"NOT_RUN_4N_FAILURE",
        "results":records}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--directory",required=True,type=pathlib.Path)
    ap.add_argument("--output",required=True,type=pathlib.Path)
    ap.add_argument("--captures",type=pathlib.Path,required=True)
    args=ap.parse_args()
    report=audit(args.directory,args.captures)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":report["status"],"methods":len(report["results"])},sort_keys=True))
if __name__=="__main__":main()
