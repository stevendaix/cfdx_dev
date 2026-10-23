"""T2 320x64 momentum-stagnation forensic diagnostic (quick level 1).

Reproduces the momentum stagnation of the finest T2 heated-channel mesh
(320x64: abort at iteration ~22, x-momentum equation residual pinned near
the rhs scale, worst cell 10561) with the solver's in-code forensic
momentum microscope enabled, and archives the dump for the audit.

Seven solves are run with the SAME case, tolerances and stop criteria as the
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
  * --debug-cell 10561 --momentum-trace-cells <see TRACE_CELLS>
    --momentum-trace-iterations 24 : the auto worst cell, plus the
    per-iteration momentum microscope trace over the inlet columns
    (i=0..4) at mid-height (j=32/33) and at the walls (j=0/63
    wall-adjacent, j=2/61 near-wall; cell id = j*320+i), bounded to
    the first 24 outer iterations. The trace also decomposes the
    x-momentum balance per face and per cell (convection, diffusion,
    pressure and body forces, boundary rhs, algebraic residual with
    the actual assembled row coefficient) and prints the effective
    relaxation factor.

Two additional A/B arms then rerun the traced solve with a controlled
momentum relaxation: alpha_u=0.4 through the diagnostic --alpha-u option,
and the solver default alpha_u=0.7 re-executed on the same binary as a
re-baseline. Mesh, boundary conditions, tolerances, gates and stop
criteria stay strictly identical; the ONLY variable is alpha_u. The
measurement is the near-wall shear-layer propagation speed: if the wall
front and the core increment scale with alpha_u per iteration, the
observed transient is a relaxation-controlled pseudo-time; if the front
arrives at the same iteration regardless of alpha_u, it follows a
physical clock independent of the relaxation.

The microscope dump (exact matrix row, Rhie-Chow face reconstruction,
pressure-correction fluxes, per-face authoritative vs reconstructed flux)
is emitted by the solver itself at the stagnation exit. The solver state,
tolerances, gates and stop criteria are untouched, and the final-flux
projection (option B) is deliberately NOT enabled here so the solve exits
through its natural stagnation path (rc=1, convergence.json/execution.json
written) instead of the projection-guard exception (rc=2) that masked the
stagnation state in run 39910523455.

Success for this driver means: the stagnation is reproduced by the five
forensic solves AND all five microscope dumps are present in the archived
stderr, and the per-iteration momentum trace of the fifth solve is
non-empty. The two alpha_u A/B arms are measurements, not gates: their
trajectories are the evidence, so an arm running past the stagnation
window or stopping earlier is a result, not an instrumentation failure.
It fails if a forensic solve converges (the symptom disappeared) or a
dump is missing (instrumentation gap). A green run is diagnostic
evidence only, never a T2 qualification verdict; the official campaign
remains scripts/t2_heated_channel_calibration.py.
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
# Traced inlet columns (i=0..4): mid-height rows j=32/33 reproduce the
# run #4 center-column trace; j=0/63 are wall-adjacent and j=2/61
# near-wall, to follow the birth of the wall shear layer (cell id =
# j*NX + i).
TRACE_CELL_ROWS = (0, 2, 32, 33, 61, 63)
TRACE_CELLS = [row * NX + i for row in TRACE_CELL_ROWS for i in range(5)]
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
    extra_args: list[str] | None = None,
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
            *(extra_args or []),
        ],
        capture_output=True,
        text=True,
        timeout=SOLVE_TIMEOUT_S,
        check=False,
    )
    diagnostics = completed.stdout + completed.stderr
    momtrace_lines = [
        line
        for line in diagnostics.splitlines()
        if line.startswith("INCOMPRESSIBLE_MOMTRACE ")
    ]
    if momtrace_lines:
        (evidence_dir / f"t2_momtrace_{NX}x{NY}_dbg{slug}_run{run_id}.log").write_text(
            "\n".join(momtrace_lines) + "\n", encoding="utf-8"
        )
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
        "momtrace_lines": len(momtrace_lines),
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
            f"  status={entry['convergence_status']} "
            f"iterations={entry['iterations']}"
        )
        lines.append(f"  stop_reason: {entry['stop_reason']}")
        lines.append(f"  json_artifacts: {', '.join(entry['json_artifacts'])}")
        lines.append(f"  momtrace_lines: {entry.get('momtrace_lines', 0)}")
    lines.append("")
    lines.append(f"alpha A/B arms: {report.get('alpha_ab', {})}")
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
        trace_args = [
            "--momentum-trace-cells",
            ",".join(str(c) for c in TRACE_CELLS),
            "--momentum-trace-iterations",
            "24",
        ]
        for slug, debug_arg, extra_args in (
            ("auto", "auto", None),
            ("cell0", "0", None),
            ("cell10241", "10241", None),
            ("cell10240", "10240", None),
            ("trace", "10561", trace_args),
            ("trace_alpha0p4", "10561", ["--alpha-u", "0.4", *trace_args]),
            ("trace_alpha0p7", "10561", trace_args),
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
                    extra_args,
                )
            )

    forensic = solves[:5]
    arms = solves[5:]
    trace_entry = forensic[-1] if forensic else {}
    if (
        all(entry["stagnation_reproduced"] for entry in forensic)
        and trace_entry.get("momtrace_lines", 0) > 0
        and all(entry.get("momtrace_lines", 0) > 0 for entry in arms)
    ):
        verdict = "stagnation_captured"
    elif any(entry["converged_marker"] for entry in forensic):
        verdict = "unexpected_convergence"
    else:
        verdict = "capture_failure"
    report = {
        "run_id": args.run_id,
        "solver": str(solver),
        "mesh": [NX, NY],
        "max_iterations": calib.MAX_ITERATIONS,
        "momentum_trace_cells": TRACE_CELLS,
        "alpha_ab": {
            "arm_alpha_u": {"trace_alpha0p4": 0.4, "trace_alpha0p7": "default"},
            "only_variable": "alpha_u",
        },
        "solves": solves,
        "verdict": verdict,
    }
    _write_evidence(out_dir, report)
    print(f"T2MOMDIAG verdict={verdict}")
    return 0 if verdict == "stagnation_captured" else 1


if __name__ == "__main__":
    sys.exit(main())
