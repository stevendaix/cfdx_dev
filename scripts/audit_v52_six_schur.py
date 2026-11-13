"""Strict same-A,b,x0 review of six *production 4N* coupled Schur attempts.

The 4N Picard saddle operator is a preconditioning reference; it is NOT the
6N Jacobian of F_phys(U,p,k,omega). This tool can never certify a BFS run.
Build the input with scripts/run_bfs_six_schur_v37.py using both
--capture-original-operator and --fingerprint-operator-v40 on <=1500 cells.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from bfs_schur_fingerprint_v40 import compare_six_fingerprints
from generate_th03_bfs_case import SCHUR_MODELS

RESIDUAL = re.compile(r"COUPLED_BLOCK_TRUE_RESIDUAL_V16[^\n]*?full_relative=([^\s]+)")
LINEAR = re.compile(r"^BFS_PHASE_V10 iteration=\d+ phase=COUPLED_LINEAR_END "
                    r"status=(\d+) iterations=(\d+) true_relative=([^\s]+) "
                    r"[^\n]*?resolved_schur=(\w+)",re.M)


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as inp:
        for chunk in iter(lambda: inp.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(campaign: dict, require_captures: bool = True) -> dict:
    names = tuple(SCHUR_MODELS)
    entries = campaign.get("results", [])
    physical = campaign.get("cases", {})
    if (len(entries) != len(names) or len(physical) != len(names) or
            {r.get("schur") for r in entries} != set(names) or
            set(physical) != set(names) or
            campaign.get("identical_physical_inputs") is not True):
        raise ValueError("six exact Schur configurations and one physical case required")
    if len({physical[name].get("physical_sha256") for name in names}) != 1:
        raise ValueError("different physical HDF5 inputs")
    if len({physical[name].get("n_cells") for name in names}) != 1:
        raise ValueError("different cell counts")
    same_fingerprint = compare_six_fingerprints(entries, names)
    hashes = []
    table = []
    for item in entries:
        method = item["schur"]
        if not item.get("operator_capture"):
            if require_captures:
                raise ValueError(f"missing exact original 4N A,b,x0 capture: {method}")
            capture = None
        else:
            capture = item["operator_capture"]
            if set(capture) != {"A", "b", "x0"}:
                raise ValueError(f"incomplete 4N capture: {method}")
            row = {}
            for name in ("A", "b", "x0"):
                if not Path(capture[name]["path"]).is_file():
                    raise ValueError(f"missing {name} capture file: {method}")
                actual = file_sha256(capture[name]["path"])
                if actual != capture[name]["sha256"]:
                    raise ValueError(f"mutated {name} capture file: {method}")
                row[name] = actual
            hashes.append(row)
        log_path = Path(item["log"])
        if not log_path.is_file() or file_sha256(log_path) != item["log_sha256"]:
            raise ValueError(f"missing or mutated execution log: {method}")
        text = log_path.read_text(encoding="utf-8", errors="replace")
        values=[]
        for raw in RESIDUAL.findall(text):
            try:
                value=float(raw)
            except ValueError as ex:
                raise ValueError(f"unparseable true residual: {method}") from ex
            if not math.isfinite(value) or value<0.0:
                raise ValueError(f"invalid true residual: {method}")
            values.append(value)
        linear=LINEAR.findall(text)
        if len(linear)>1:
            raise ValueError(f"multiple frozen 4N Krylov observations: {method}")
        if linear:
            status, iterations, relative, resolved=linear[0]
            if resolved!=method:
                raise ValueError(f"silently substituted Schur method: {method} => {resolved}")
            try: relative_value=float(relative)
            except ValueError as ex:
                raise ValueError(f"malformed Krylov true residual: {method}") from ex
            if not math.isfinite(relative_value) or relative_value<0.0:
                raise ValueError(f"nonfinite Krylov true residual: {method}")
            true_relative=relative_value
            solver_iterations=int(iterations)
            solver_converged=(int(status)==0)
        else:
            true_relative=values[-1] if values else None
            solver_iterations=None
            solver_converged=False
        table.append({"schur":method,"execution":item.get("execution"),
                      "first_blocker":item.get("first_blocker",{}).get("first_blocker"),
                      "wall_seconds":item.get("wall_seconds"),
                      "last_4n_true_relative":true_relative,
                      "4n_linear_iterations":solver_iterations,
                      "4n_linear_solver_converged":solver_converged,
                      "nonlinear_converged":item.get("nonlinear_convergence_proven") is True,
                      "6n_jfnk_preconditioner_qualified":False})
    same_sha=(len(hashes)==6 and
              all(len({row[name] for row in hashes})==1
                  for name in ("A","b","x0")))
    return {"schema":"CFDX-V52-SIX-SCHUR-COMMON-4N-AUDIT/1",
            "n_cells":next(iter(physical.values()))["n_cells"],
            "same_physical_HDF5":True,
            "all_six_same_original_4N_A_b_x0_SHA256":same_sha,
            "all_six_same_original_4N_A_b_x0_FNV64":
                same_fingerprint["all_six_same_A_b_x0_fingerprint"],
            "matrices_are_6N_Jacobian":False,
            "globalized_6N_preconditioners_compared":False,
            "physically_qualified":False,
            "status":"COMPARABLE_4N_ONLY" if same_sha and
                same_fingerprint["all_six_same_A_b_x0_fingerprint"] else
                "INCOMPARABLE_OR_INCOMPLETE",
            "methods":sorted(table,key=lambda row: names.index(row["schur"]))}


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    try:
        result=audit(json.loads(args.campaign.read_text(encoding="utf-8")))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",
                               encoding="utf-8")
        print("V52 six Schur: "+result["status"]+"; 6N JFNK=NOT_QUALIFIED")
        return 0 if result["status"]=="COMPARABLE_4N_ONLY" else 2
    except (OSError,ValueError,KeyError,TypeError,json.JSONDecodeError) as error:
        print(f"V52 six Schur proof invalid: {error}",file=sys.stderr)
        return 3

if __name__=="__main__":
    raise SystemExit(main())
