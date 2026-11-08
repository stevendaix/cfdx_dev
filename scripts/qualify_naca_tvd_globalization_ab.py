"""Fail-closed, hash-traceable NACA SAME-MESH TVD Picard experiment driver.

An iteration-convergence certificate is NOT a physical Cd/Cf/GCI certificate.
Only comparable cases (same binary/mesh/relaxations/linear tolerances) are run.
Never downgrades existing strict nonlinear or face-residual solver gates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

CASES: dict[str, dict[str, str]] = {
    "baseline": {},
    "memory_raw_025": {"CFDX_NACA_TVD_DEFECT_MEMORY": "0.25"},
    "memory_recursive_025": {
        "CFDX_NACA_TVD_DEFECT_MEMORY": "0.25",
        "CFDX_NACA_TVD_MEMORY_RECURSIVE": "1",
    },
    "memory_oscillation_025": {
        "CFDX_NACA_TVD_DEFECT_MEMORY": "0.25",
        "CFDX_NACA_TVD_MEMORY_OSCILLATION": "1",
    },
    "secant_050": {"CFDX_NACA_TVD_IMPLICIT_SECANT": "0.5"},
    "globalization_nonmonotone": {
        "CFDX_NACA_PICARD_GLOBALIZATION": "1",
        "CFDX_NACA_PICARD_WINDOW": "8",
        "CFDX_NACA_PICARD_GROWTH": "1.05",
    },
}
GATES = {
    "momentum_residual_rhs_relative": 1e-8,
    "momentum_active_component_residual_max": 1e-8,
    "continuity_normalized": 1e-8,
    "pressure_relative_residual": 1e-8,
    "velocity_change_inf": 1e-8,
    "pressure_change_inf": 1e-8,
}
ITER_LINE = "INCOMPRESSIBLE_ITER "
MEM_LINE = "TVD_DEFECT_MEMORY "


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def fields(line: str) -> dict[str, str]:
    return dict(re.findall(r"([A-Za-z_][A-Za-z_0-9]*)=([^\s]+)", line))


def diagnose(stdout: str, stderr: str, exit_code: int,
             mode: str, timed_out: bool = False) -> dict[str, Any]:
    """Parse an actual NACA executable output; never infer PASS
       from a low Krylov residual alone."""
    traces=[fields(x) for x in stderr.splitlines() if ITER_LINE in x]
    if not traces:
        return {"solver_contract_pass":False, "reason":"MISSING_PHYSICAL_ITERATION",
                "exit_code":exit_code,"iterations_observed":0}
    last=traces[-1]
    observed: dict[str, float] = {}
    for k in GATES:
        try:
            observed[k]=float(last[k])
        except (ValueError,KeyError):
            return {"solver_contract_pass":False,"reason":f"MISSING_METRIC_{k}",
                    "exit_code":exit_code,"iterations_observed":len(traces)}
        if not (observed[k] >= 0 and observed[k] < float("inf")):
            return {"solver_contract_pass":False,"reason":f"INVALID_METRIC_{k}",
                    "exit_code":exit_code,"iterations_observed":len(traces)}
    if mode.startswith("memory_"):
        memory=[fields(x) for x in stderr.splitlines() if MEM_LINE in x]
        if not memory or "source_relative" not in memory[-1]:
            return {"solver_contract_pass":False,"reason":"MISSING_TVD_SOURCE_GATE",
                    "exit_code":exit_code,"metrics":observed}
        try:
            observed["tvd_source_consistency_relative"]=float(
                memory[-1]["source_relative"])
        except ValueError:
            return {"solver_contract_pass":False,"reason":"INVALID_TVD_SOURCE_GATE",
                    "exit_code":exit_code,"metrics":observed}
    passed_labels = ("NACA0012_SINGLE: PASS" in stdout and
                     "NACA0012_QUALIFICATION: FAIL" not in stderr and
                     exit_code == 0 and not timed_out)
    gates_ok = (all(observed[k]<=v for k,v in GATES.items()) and
                last.get("high_order_fraction") in ("1","1.0") and
                last.get("full_order_eligible")=="1" and
                observed.get("tvd_source_consistency_relative",0.0) <= 1e-8)
    passed=bool(passed_labels and gates_ok)
    return {"solver_contract_pass":passed,
            "physical_reference_qualified":False,
            "reason": ("ITERATIVE_CONVERGENCE_ONLY" if passed else
                       "INCOMPLETE_OR_NONCONVERGED"),
            "exit_code":exit_code,"timed_out":timed_out,
            "last_iteration":int(last.get("iter","0")),
            "iterations_observed":len(traces),"metrics":observed,
            "full_order_fraction":last.get("high_order_fraction"),
            "full_order_eligible":last.get("full_order_eligible")}


def run_campaign(*, binary: Path, mesh: Path, output: Path,
                 alpha_u: float, alpha_p: float, iterations: int,
                 seconds: float, case_names: list[str]) -> dict[str, Any]:
    binary=binary.resolve(strict=True)
    mesh=mesh.resolve(strict=True)
    output=output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    result:dict[str,Any]={
        "schema":"cfdx_naca_picard_ab.v1",
        "scope":"NACA N64 same-mesh numerical investigation",
        "criteria":GATES, "physical_reference_qualified":False,
        "binary_sha256":sha256_file(binary), "mesh_sha256":sha256_file(mesh),
        "alpha_u":alpha_u,"alpha_p":alpha_p,"max_iterations":iterations,"cases":{},
    }
    for case_name in case_names:
        cfg=CASES[case_name]
        location=output/case_name
        location.mkdir(parents=True,exist_ok=True)
        clean={k:v for k,v in os.environ.items() if not k.startswith("CFDX_NACA_")}
        clean["CFDX_NACA_CONVECTION_TVD"]="1"
        clean.update(cfg)
        cmd=[str(binary),"--single",str(mesh),str(alpha_u),str(alpha_p),str(iterations)]
        timed_out=False
        try:
            proc=subprocess.run(cmd,cwd=location,env=clean,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                text=True,encoding="utf-8",errors="replace",
                                timeout=seconds,check=False)
            stdout,stderr,return_code=proc.stdout,proc.stderr,proc.returncode
        except subprocess.TimeoutExpired as error:
            timed_out=True
            stdout=error.stdout or b""
            stderr=error.stderr or b""
            if isinstance(stdout,bytes):stdout=stdout.decode("utf-8","replace")
            if isinstance(stderr,bytes):stderr=stderr.decode("utf-8","replace")
            return_code=124
        (location/"stdout.log").write_text(stdout,encoding="utf-8")
        (location/"stderr.log").write_text(stderr,encoding="utf-8")
        result["cases"][case_name]={"env":{"CFDX_NACA_CONVECTION_TVD":"1",**cfg},
            **diagnose(stdout,stderr,return_code,case_name,timed_out),
            "stdout_sha256":sha256_file(location/"stdout.log"),
            "stderr_sha256":sha256_file(location/"stderr.log")}
        (output/"summary.json").write_text(
            json.dumps(result,indent=2,sort_keys=True)+"\n")
    result["all_solver_contracts_pass"]=all(
        c["solver_contract_pass"] for c in result["cases"].values())
    (output/"summary.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    return result


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--mesh",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--alpha-u",type=float,default=0.6)
    parser.add_argument("--alpha-p",type=float,default=0.25)
    parser.add_argument("--iterations",type=int,default=200)
    parser.add_argument("--timeout",type=float,default=90.0)
    parser.add_argument("--case",dest="cases",choices=tuple(CASES),action="append")
    args=parser.parse_args()
    if (not (0<args.alpha_u<=1 and 0<args.alpha_p<=1) or
            args.iterations<1 or not (0<args.timeout<1e5)):
        parser.error("invalid relaxations, iteration budget or runtime")
    output=run_campaign(binary=args.binary,mesh=args.mesh,output=args.output,
                        alpha_u=args.alpha_u,alpha_p=args.alpha_p,
                        iterations=args.iterations,seconds=args.timeout,
                        case_names=args.cases or
                        ["baseline","memory_raw_025","secant_050"])
    for name,entry in output["cases"].items():
        print(name,entry["reason"],entry.get("metrics",{}).get(
            "momentum_active_component_residual_max"))
    print("All numerical gates pass:",output["all_solver_contracts_pass"])
    # A diagnostic campaign with no validated solution must fail in CI.
    return 0 if output["all_solver_contracts_pass"] else 2

if __name__ == "__main__":
    sys.exit(main())

