"""T2 320x64 momentum-stagnation forensic diagnostic (quick level 1).

Reproduces the momentum stagnation of the finest T2 heated-channel mesh
(320x64: abort at iteration ~22, x-momentum equation residual pinned near
the rhs scale, worst cell 10561) with the solver's in-code forensic
momentum microscope enabled, and archives the dump for the audit.

Four solves are run with the SAME case, tolerances and stop criteria as the
calibration campaign (module t2_heated_channel_calibration is imported,
never duplicated):

  * --debug-cell auto : the independently reassembled worst momentum
    residual cell (10561 in runs #8-#10),
  * --debug-cell 0 : the pressure-gauge cell at the inlet/bottom corner,
    the worst mass-continuity cell since iteration 3.
  * --debug-cell 10241 : cell (i=1, j=32), one column past the inlet at
    mid-height, a member of the worst-cell cluster of run 39939509005,
  * --debug-cell 10240 : cell (i=0, j=32), the inlet column at
    mid-height, adjacent to the boundary faces that dominate the
    boundary-driven share of the pinned momentum residual.

The microscope dump (exact matrix row, Rhie-Chow face reconstruction,
pressure-correction fluxes, per-face authoritative vs reconstructed flux)
is emitted by the solver itself at the stagnation exit. The solver state,
tolerances, gates and stop criteria are untouched, and the final-flux
projection (option B) is deliberately NOT enabled here so the solve exits
through its natural stagnation path (rc=1, convergence.json/execution.json
written) instead of the projection-guard exception (rc=2) that masked the
stagnation state in run 39910523455.

Success for this driver means: the stagnation is reproduced AND all four
microscope dumps are present in the archived stderr. It fails if the
solve converges (the symptom disappeared) or a dump is missing
(instrumentation gap). A green run is diagnostic evidence only, never a
T2 qualification verdict; the official campaign remains
scripts/t2_heated_channel_calibration.py.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

NX, NY = 320, 64
SOLVE_TIMEOUT_S = 3600
MICROSCOPE_BEGIN = "=== CFDX MOMENTUM MICROSCOPE cell="
MICROSCOPE_END = "=== END CFDX MOMENTUM MICROSCOPE ==="
STAGNATION_REASON_SNIPPET = "stagnation window"


def _load_calibration_module():
    """Import the sibling calibration driver as a module (no sys.path)."""
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(
        "t2_heated_channel_calibration", here / "t2_heated_channel_calibration.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load t2_heated_channel_calibration")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _scalar(artifact, *keys):
    """First present key of a json artifact; None-safe for absent files."""
    if not isinstance(artifact, dict):
        return None
    for key in keys:
        if artifact.get(key) is not None:
            return artifact[key]
    return None


def _run_solve(
    calib,
    solver: Path,
    slug: str,
    debug_arg: str,
    run_dir: Path,
    evidence_dir: Path,
    run_id: str,
) -> dict:
    """One 320x64 forensic solve; never raises, evidence always archived."""
    run_dir.mkdir(parents=True, exist_ok=True)
    case_path = run_dir / f"t2_heated_channel_{NX}x{NY}.cfdx.h5"
    calib._write_case(case_path, NX, NY)
    completed = subprocess.run(
        [
            str(solver),
            "--mesh",
            str(case_path),
            "--output-dir",
            str(run_dir),
            "--iterations",
            str(calib.MAX_ITERATIONS),
            "--debug-cell",
            debug_arg,
        ],
        capture_output=True,
        text=True,
        timeout=SOLVE_TIMEOUT_S,
        check=False,
    )
    diagnostics = completed.stdout + completed.stderr
    (run_dir / f"solver_{slug}.log").write_text(diagnostics, encoding="utf-8")
    (evidence_dir / f"t2_solver_tail_{NX}x{NY}_dbg{slug}_run{run_id}.log").write_text(
        diagnostics[-200000:], encoding="utf-8"
    )
    (evidence_dir / f"t2_solver_stdout_{NX}x{NY}_dbg{slug}_run{run_id}.log").write_text(
        completed.stdout[-400000:], encoding="utf-8"
    )
    artifacts: dict = {}
    for artifact in sorted(run_dir.glob("*.json")):
        text = artifact.read_text(encoding="utf-8")
        target = evidence_dir / (
            f"t2_{artifact.stem}_{NX}x{NY}_dbg{slug}_run{run_id}.json"
        )
        target.write_text(text, encoding="utf-8")
        try:
            artifacts[artifact.stem] = json.loads(text)
        except ValueError:
            artifacts[artifact.stem] = None

    convergence = artifacts.get("convergence")
    execution = artifacts.get("execution")
    status = _scalar(execution, "convergence_status", "status")
    if status is None:
        status = _scalar(convergence, "convergence_status", "status")
    iterations = _scalar(execution, "iterations")
    if iterations is None:
        iterations = _scalar(convergence, "iterations")
    reason_value = _scalar(convergence, "convergence_reason", "reason")
    if reason_value is None:
        reason_value = _scalar(execution, "convergence_reason", "reason")
    reason = str(reason_value) if reason_value is not None else ""
    match = re.search(re.escape(MICROSCOPE_BEGIN) + r"(\d+)", diagnostics)
    entry: dict = {
        "debug_cell_arg": debug_arg,
        "returncode": completed.returncode,
        "converged_marker": "Converged YES" in diagnostics,
        "microscope_present": MICROSCOPE_BEGIN in diagnostics,
        "microscope_complete": MICROSCOPE_END in diagnostics,
        "microscope_cell": int(match.group(1)) if match else None,
        "convergence_status": status,
        "iterations": iterations,
        "stop_reason": reason,
        "json_artifacts": sorted(artifacts),
    }
    entry["stagnation_reproduced"] = (
        completed.returncode == 1
        and "Converged NO" in diagnostics
        and MICROSCOPE_BEGIN in diagnostics
        and STAGNATION_REASON_SNIPPET in reason
    )
    print(
        f"T2MOMDIAG_SOLVE slug={slug} debug_cell={debug_arg} "
        f"rc={completed.returncode} microscope={entry['microscope_present']} "
        f"cell={entry['microscope_cell']} "
        f"stagnation={entry['stagnation_reproduced']}"
    )
    return entry


def _write_evidence(out_dir: Path, report: dict) -> None:
    run_id = report["run_id"]
    json_path = out_dir / f"t2_momdiag_{run_id}.json"
    txt_path = out_dir / f"t2_momdiag_{run_id}.txt"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    lines: list[str] = []
    lines.append("T2 320x64 momentum-stagnation forensic diagnostic (level 1)")
    lines.append(f"run_id: {run_id}")
    lines.append(f"solver: {report['solver']}")
    lines.append("mesh: 320x64, max_iterations=3000 (calibration module reused)")
    lines.append("no gate, tolerance, or stop criterion modified; no flux projection")
    lines.append("")
    for entry in report["solves"]:
        lines.append(
            f"solve --debug-cell {entry['debug_cell_arg']}: "
            f"rc={entry['returncode']} "
            f"microscope_present={entry['microscope_present']} "
            f"microscope_complete={entry['microscope_complete']} "
            f"microscope_cell={entry['microscope_cell']} "
            f"stagnation_reproduced={entry['stagnation_reproduced']}"
        )
        lines.append(
            f"  status={entry['convergence_status']} iterations={entry['iterations']}"
        )
        lines.append(f"  stop_reason: {entry['stop_reason']}")
        lines.append(f"  json_artifacts: {', '.join(entry['json_artifacts'])}")
    lines.append("")
    lines.append(f"verdict: {report['verdict']}")
    lines.append(
        "interpretation: a green run captures forensic evidence only; it"
        " is never a T2 qualification verdict."
    )
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    calib = _load_calibration_module()
    calib._ensure_cfdx_package()
    parser = argparse.ArgumentParser(
        description="T2 320x64 momentum-stagnation forensic diagnostic"
    )
    parser.add_argument("--solver", required=True, help="production solver binary")
    parser.add_argument("--output-dir", required=True, help="evidence directory")
    parser.add_argument("--run-id", default="local", help="diagnostic run id")
    args = parser.parse_args()

    solver = Path(args.solver)
    if not solver.is_file():
        print(f"T2MOMDIAG FAIL solver_missing={solver}")
        return 1
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    solves: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="t2_momdiag_") as tmp:
        work = Path(tmp)
        for slug, debug_arg in (
            ("auto", "auto"),
            ("cell0", "0"),
            ("cell10241", "10241"),
            ("cell10240", "10240"),
        ):
            solves.append(
                _run_solve(
                    calib,
                    solver,
                    slug,
                    debug_arg,
                    work / f"solve_{slug}",
                    out_dir,
                    args.run_id,
                )
            )

    if all(entry["stagnation_reproduced"] for entry in solves):
        verdict = "stagnation_captured"
    elif any(entry["converged_marker"] for entry in solves):
        verdict = "unexpected_convergence"
    else:
        verdict = "capture_failure"
    report = {
        "run_id": args.run_id,
        "solver": str(solver),
        "mesh": [NX, NY],
        "max_iterations": calib.MAX_ITERATIONS,
        "solves": solves,
        "verdict": verdict,
    }
    _write_evidence(out_dir, report)
    print(f"T2MOMDIAG verdict={verdict}")
    return 0 if verdict == "stagnation_captured" else 1


if __name__ == "__main__":
    sys.exit(main())
