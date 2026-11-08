"""Execute a physically gated one-level NACA study; keep logs if it fails."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path


def run_case(mesh: Path, output: Path, exe: Path, maxit: int, seconds: int,
             initial_state: Path | None = None,
             export_converged_state: Path | None = None,
             initial_face_flux: Path | None = None,
             export_converged_face_flux: Path | None = None):
    if initial_face_flux is not None and initial_state is None:
        raise ValueError(
            "face flux import requires a paired same-mesh cell state")
    if initial_face_flux is not None and not initial_face_flux.is_file():
        raise FileNotFoundError(
            f"not a valid face-flux restart file: {initial_face_flux}")
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({
        "CFDX_NACA_CONVECTION_VANLEER": "1",
        "CFDX_NACA_LINEAR_TOLERANCE": "1e-12",
        "CFDX_NACA_WALL_TRACTION_FACE": "1",
        "CFDX_NACA_WALL_GRADIENT_WLS": "1",
        "CFDX_NACA_CONVERGED_FORCE_CSV": str(output / "force.csv"),
        "CFDX_NACA_PHYSICAL_HISTORY_CSV": str(output / "history.csv"),
        "CFDX_NACA_FACE_PROFILES_CSV": str(output / "profiles.csv"),
        "CFDX_NACA_CHECKPOINT_STATE_CSV": str(output / "checkpoint.csv"),
        "CFDX_NACA_CHECKPOINT_EVERY": "300",
    })
    if initial_state is not None:
        if not initial_state.is_file():
            raise FileNotFoundError(
                f"not a valid restart initial state: {initial_state}")
        env["CFDX_NACA_IMPORT_INITIAL_STATE_CSV"] = str(initial_state.resolve())
    if export_converged_state is not None:
        env["CFDX_NACA_EXPORT_CONVERGED_STATE_CSV"] = str(
            export_converged_state.resolve())
    if initial_face_flux is not None:
        env["CFDX_NACA_IMPORT_FACE_FLUX_CSV"] = str(initial_face_flux.resolve())
    if export_converged_face_flux is not None:
        env["CFDX_NACA_CONVERGED_FACE_FLUX_CSV"] = str(
            export_converged_face_flux.resolve())
    start = time.monotonic()
    with (output / "run.log").open("w") as log:
        proc = subprocess.run(
            ["timeout", f"{seconds}s", str(exe), "--single", str(mesh),
             "0.6", "0.25", str(maxit)],
            env=env, stdout=log, stderr=subprocess.STDOUT,
            check=False,
        )
    elapsed = time.monotonic() - start
    (output / "run_summary.json").write_text(
        json.dumps({
            "mesh": str(mesh), "scheme": "TVD_VANLEER",
            "linear_tolerance": 1e-12, "wall_traction": "FACE_WLS",
            "nonlinear_tolerance": 1e-8, "max_iter": maxit,
            "runtime_seconds": round(elapsed, 2),
            "exit_code": proc.returncode,
            "certificate_present": (output / "force.csv").exists(),
            "history_present": (output / "history.csv").exists(),
            "restart_state_used_as_initial_guess": (
                str(initial_state.resolve()) if initial_state else None),
            "restart_face_flux_input": (
                str(initial_face_flux.resolve()) if initial_face_flux else None),
            "restart_fidelity": (
                "PAIRED_CONSERVATIVE_FACE_FLUX" if initial_face_flux
                else "CELL_STATE_ONLY_FACE_FLUX_RECONSTRUCTED" if initial_state
                else "COLD_START"),
            "converged_face_flux_present": bool(
                export_converged_face_flux
                and export_converged_face_flux.is_file()),
        }, indent=2) + "\n")
    print(f"{mesh.stem} exit={proc.returncode} runtime={elapsed:.2f}s "
          f"certificate={(output / 'force.csv').exists()}", flush=True)
    return proc.returncode


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--exe", type=Path,
        default=Path("/mnt/data/cfdx_work/build/test_naca0012_qualification"))
    parser.add_argument("--maxit", type=int, default=2600)
    parser.add_argument("--seconds", type=int, default=1100)
    parser.add_argument("--initial-state", type=Path)
    parser.add_argument("--export-converged-state", type=Path)
    parser.add_argument(
        "--initial-face-flux", type=Path,
        help="same-mesh conservative face flux paired with --initial-state")
    parser.add_argument(
        "--export-converged-face-flux", type=Path,
        help="write Rhie-Chow face flux only after physical convergence")
    args = parser.parse_args()
    raise SystemExit(run_case(
        args.mesh, args.output_dir, args.exe, args.maxit, args.seconds,
        args.initial_state, args.export_converged_state,
        args.initial_face_flux, args.export_converged_face_flux))
