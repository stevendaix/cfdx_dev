"""Execute real CFDX TH benchmark cases, never synthetic solver profiles.

The supplied campaign JSON must point at actual .cfdx.h5 cases. Unsupported
physics (LES/DNS, periodic channel, SST+buoyancy/energy, transient statistics)
is *blocked before execution*. This launcher does not grant physical validity;
users must separately provide independent, pinned reference profiles and run
qualify_thermohydraulic_benchmarks.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

try:
    from scripts.qualify_thermohydraulic_benchmarks import SCENARIOS, MAX_MASS_IMBALANCE, MAX_TRUE_RESIDUAL
except ImportError:
    from qualify_thermohydraulic_benchmarks import SCENARIOS, MAX_MASS_IMBALANCE, MAX_TRUE_RESIDUAL


class CampaignError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CampaignError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def positive_number(value: Any, name: str) -> float:
    require(isinstance(value, (float, int)) and not isinstance(value, bool),
            f"{name}: expected a positive number")
    result = float(value)
    require(math.isfinite(result) and result > 0, f"{name}: must be positive and finite")
    return result


def resolve_relative(root: Path, value: str, name: str) -> Path:
    require(isinstance(value, str) and bool(value), f"{name}: empty path")
    file = (root / value).resolve()
    require(file.is_relative_to(root.resolve()), f"{name}: escapes campaign directory")
    return file


def preflight(campaign: dict[str, Any], root: Path, case_filter: set[str] | None) -> dict[str, Any]:
    require(campaign.get("schema_version") == 1, "campaign schema_version must be 1")
    cases = campaign.get("cases")
    require(isinstance(cases, list) and cases, "campaign cases must be nonempty")
    ids = [item.get("id") for item in cases if isinstance(item, dict)]
    require(len(ids) == len(cases) and len(set(ids)) == len(ids),
            "duplicate or malformed benchmark case ID")
    require(set(ids).issubset(SCENARIOS), "unknown TH benchmark case ID")
    selected = [item for item in cases if case_filter is None or item["id"] in case_filter]
    require(selected, "no selected TH benchmark case in campaign")
    report: dict[str, Any] = {"schema_version": 1, "phase": "capability-and-input-preflight",
                              "never_a_physical_validation": True, "cases": []}
    for case in selected:
        case_id = case["id"]
        model = case.get("model")
        require(model in ("RANS_SST", "RSM", "LES", "DNS"),
                f"{case_id}: model must be declared explicitly")
        scenarios = case.get("scenarios")
        require(isinstance(scenarios, list), f"{case_id}: scenarios absent")
        scenario_ids = [s.get("id") for s in scenarios if isinstance(s, dict)]
        require(set(scenario_ids) == set(SCENARIOS[case_id]) and
                len(scenario_ids) == len(scenarios),
                f"{case_id}: missing or duplicate benchmark scenarios")
        blockers: list[str] = []
        # These are verified against the production executable's current
        # supported paths. They cannot be bypassed by toggling a JSON flag.
        if case_id == "TH-01":
            blockers.append("periodic 3D, mass-flux-controlled channel DNS/LES is absent")
        if case_id == "TH-07":
            blockers.append("SST+energy+buoyancy production coupling is absent")
        if case_id == "TH-08":
            blockers.append("transient SST/LES energy and sampling/spectra are absent")
        if model in ("LES", "DNS", "RSM"):
            blockers.append(f"{model} transport is not available in the production driver")
        if model == "RANS_SST":
            blockers.append("SST remains experimental and not physically qualified")
        ready_inputs = True
        for scenario in scenarios:
            sid = scenario["id"]
            runs = scenario.get("runs")
            require(isinstance(runs, list) and len(runs) >= 3,
                    f"{case_id}/{sid}: three mesh levels are mandatory")
            ids_seen: set[str] = set()
            last_cells, last_h = 0.0, math.inf
            for i, run in enumerate(runs):
                label = f"{case_id}/{sid}/mesh[{i}]"
                require(isinstance(run, dict), f"{label}: malformed entry")
                mid = run.get("mesh_id")
                require(isinstance(mid, str) and re.fullmatch(r"[a-zA-Z0-9_-]+",mid) is not None,
                        f"{label}: mesh_id must be filesystem-safe")
                require(mid not in ids_seen, f"{label}: duplicate mesh_id")
                ids_seen.add(mid)
                cells = positive_number(run.get("cells"),f"{label}.cells")
                h = positive_number(run.get("characteristic_h"),f"{label}.characteristic_h")
                require(cells > last_cells and h < last_h and cells >= 100,
                        f"{label}: mesh refinement ordering invalid")
                last_cells, last_h = cells, h
                iters = run.get("max_iterations")
                require(type(iters) is int and 0 < iters <= 10_000_000,
                        f"{label}: max_iterations out of range")
                case_file = resolve_relative(root,run.get("case_file"),f"{label}.case_file")
                require(case_file.name.endswith(".cfdx.h5"),
                        f"{label}: production benchmarks require a .cfdx.h5 case")
                if not case_file.is_file():
                    ready_inputs = False
                    blockers.append(f"missing real .cfdx.h5 input {case_file}")
        # All RANS_SST inputs remain eligible to *execute* only TH-03 as an
        # explicitly experimental numerical investigation, never qualification.
        runnable = ready_inputs and case_id == "TH-03" and model == "RANS_SST"
        report["cases"].append({"id":case_id, "model":model,
            "inputs_present":ready_inputs, "executable_experimental":runnable,
            "physical_qualification":"NOT_QUALIFIED", "blockers":list(dict.fromkeys(blockers))})
    report["status"] = ("EXPERIMENTAL_RUN_INPUTS_READY" if any(
        c["executable_experimental"] for c in report["cases"]) else "BLOCKED")
    return report


def run_one(solver: Path, campaign_root: Path, out: Path,
            case_id: str, scenario: dict[str, Any], run: dict[str, Any],
            timeout: int) -> dict[str, Any]:
    sid, mid = scenario["id"], run["mesh_id"]
    run_dir = out / case_id / sid / mid
    run_dir.mkdir(parents=True, exist_ok=True)
    case_file = resolve_relative(campaign_root,run["case_file"],"case_file")
    args = [str(solver),"--mesh",str(case_file),"--output-dir",str(run_dir),
            "--iterations",str(run["max_iterations"]),
            "--experimental-sst","--benchmark-export"]
    # Deliberately no arbitrary user-controlled argv or shell invocation.
    log = run_dir / "solver.log"
    try:
        with log.open("wb") as stream:
            proc = subprocess.run(args, stdout=stream,stderr=subprocess.STDOUT,
                                  timeout=timeout, check=False)
        exit_code = proc.returncode
        error = "" if exit_code == 0 else f"CFDX process exited with code {exit_code}"
    except subprocess.TimeoutExpired:
        exit_code, error = -1, f"CFDX run exceeded the {timeout}s budget"
    execution_file = run_dir / "execution.json"
    convergence_file = run_dir / "convergence.json"
    flux_file = run_dir / "benchmark_observables" / "mass_conservation.json"
    actual = {}
    if exit_code == 0:
        try:
            execution = json.loads(execution_file.read_text(encoding="utf-8"))
            convergence = json.loads(convergence_file.read_text(encoding="utf-8"))
            conservation = json.loads(flux_file.read_text(encoding="utf-8"))
            require(execution.get("format") == "CFDX-EXECUTION" and
                    execution.get("process_exit_code") == 0 and execution.get("converged") is True,
                    "execution.json contradicts process exit")
            require(convergence.get("format") == "CFDX-CONVERGENCE" and
                    convergence.get("converged") is True and convergence.get("history"),
                    "convergence.json is incomplete")
            require(conservation.get("format") == "CFDX-THERMOHYDRAULIC-CONSERVATION",
                    "mass_conservation.json missing")
            last = convergence["history"][-1]
            momentum = float(last["momentum_equation_residual_rhs_relative"])
            continuity = float(last["continuity_normalized"])
            require(math.isfinite(momentum) and momentum >= 0 and
                    math.isfinite(continuity) and continuity >= 0,
                    "invalid nonlinear residual metrics")
            residual = max(momentum, continuity)
            mass = float(conservation["global_mass_imbalance_relative"])
            local_mass = float(conservation["maximum_cell_mass_imbalance_relative"])
            require(math.isfinite(mass) and mass >= 0 and
                    math.isfinite(local_mass) and local_mass >= 0,
                    "invalid mass conservation measures")
            require(residual <= MAX_TRUE_RESIDUAL and
                    mass <= MAX_MASS_IMBALANCE and local_mass <= MAX_MASS_IMBALANCE,
                    "native residual, global or local mass balance exceeds fixed gate")
            require(conservation.get("cells") == int(run["cells"]),
                    "native CFDX cell count disagrees with declared input mesh")
            actual = {"iterations":execution["iterations"],
                      "true_residual_max":residual,
                      "mass_imbalance_relative":mass,
                      "maximum_cell_mass_imbalance_relative":local_mass,
                      "files_sha256":{
                         str(path.relative_to(run_dir)):sha256_file(path)
                         for path in (execution_file,convergence_file,flux_file,
                                      run_dir / "benchmark_observables" / "cell_fields.csv",
                                      run_dir / "benchmark_observables" / "face_fluxes.csv")}}
        except (OSError,KeyError,ValueError,TypeError,CampaignError,OverflowError) as exc:
            error = f"invalid solver-owned evidence: {exc}"
    result: dict[str, Any] = {"case_id":case_id,"scenario_id":sid,"mesh_id":mid,
        "case_file_sha256":sha256_file(case_file),
        "solver_binary_sha256":sha256_file(solver),"solver_log_sha256":sha256_file(log),
        "exit_code":exit_code,"solver_executed":True,
        "physical_qualification":"NOT_QUALIFIED",
        "status":"EXECUTED_UNQUALIFIED" if not error else "INCOMPLETE_OR_FAIL",
        "error":error, **actual}
    (run_dir / "run_evidence.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign",required=True,type=Path)
    parser.add_argument("--case",action="append",choices=tuple(SCENARIOS))
    parser.add_argument("--check-only",action="store_true")
    parser.add_argument("--solver",type=Path)
    parser.add_argument("--output",type=Path,default=Path("build/thermohydraulic/native-runs"))
    parser.add_argument("--timeout",type=int,default=7200)
    args = parser.parse_args(argv)
    result_dir = args.output.resolve()
    result_dir.mkdir(parents=True,exist_ok=True)
    report_path = result_dir / "campaign-report.json"
    try:
        root = args.campaign.resolve().parent
        campaign = json.loads(args.campaign.read_text(encoding="utf-8"))
        pre = preflight(campaign,root,set(args.case) if args.case else None)
        report: dict[str,Any] = {"preflight":pre,"runs":[],
            "never_a_physical_validation":True,"status":"BLOCKED"}
        if not args.check_only:
            require(args.solver is not None and args.solver.resolve().is_file(),
                    "--solver must point to an existing cfdx_production_solver binary")
            require(type(args.timeout) is int and args.timeout > 0,"timeout must be positive")
            selected = [c for c in campaign["cases"] if not args.case or c["id"] in args.case]
            for item in selected:
                eligibility = next(v for v in pre["cases"] if v["id"] == item["id"])
                if not eligibility["executable_experimental"]:
                    continue
                for scenario in item["scenarios"]:
                    for run in scenario["runs"]:
                        report["runs"].append(run_one(
                            args.solver.resolve(),root,result_dir,item["id"],scenario,run,args.timeout))
            report["status"] = ("EXECUTED_UNQUALIFIED" if report["runs"] and
                all(r["status"] == "EXECUTED_UNQUALIFIED" for r in report["runs"])
                else "INCOMPLETE_OR_FAIL")
        else:
            report["status"] = pre["status"]
        report_path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        print(json.dumps(report,indent=2,allow_nan=False))
        # Check-only is allowed to exit zero after reporting blockers; the
        # physical execution path never exits zero for incomplete scenarios.
        return 0 if args.check_only else (0 if report["status"] == "EXECUTED_UNQUALIFIED" else 2)
    except (CampaignError,OSError,json.JSONDecodeError,ValueError,TypeError) as exc:
        report = {"status":"INCOMPLETE_OR_FAIL","error":str(exc),
                  "physical_qualification":"NOT_QUALIFIED"}
        report_path.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print(json.dumps(report,indent=2),file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
