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
             export_converged_face_flux: Path | None = None,
             initial_paired_checkpoint: Path | None = None,
             checkpoint_every: int = 300,
             stagnation_diagnostic_only: bool = False):
    if checkpoint_every < 10 or checkpoint_every > 10000:
        raise ValueError("checkpoint interval must be in [10,10000]")
    if initial_paired_checkpoint is not None and (
            initial_face_flux is not None or initial_state is not None):
        raise ValueError(
            "paired manifest restart is exclusive with legacy state/flux inputs")
    if (initial_paired_checkpoint is not None
            and not initial_paired_checkpoint.is_file()):
        raise FileNotFoundError(
            f"paired restart manifest is not readable: {initial_paired_checkpoint}")
    if initial_face_flux is not None and initial_state is None:
        raise ValueError(
            "face flux import requires a paired same-mesh cell state")
    if initial_face_flux is not None and not initial_face_flux.is_file():
        raise FileNotFoundError(
            f"not a valid face-flux restart file: {initial_face_flux}")
    output.mkdir(parents=True, exist_ok=True)
    # Same directory cannot simultaneously be the source and destination.
    manifest = output / "checkpoint.manifest"
    if (initial_paired_checkpoint is not None
            and manifest.resolve() == initial_paired_checkpoint.resolve()):
        raise ValueError("cannot overwrite the active restart input manifest")
    # A stale checkpoint must never be reported as a new checkpoint after
    # failure.
    manifest.unlink(missing_ok=True)
    env = os.environ.copy()
    env.update({
        "CFDX_NACA_CONVECTION_VANLEER": "1",
        "CFDX_NACA_LINEAR_TOLERANCE": "1e-12",
        "CFDX_NACA_WALL_TRACTION_FACE": "1",
        "CFDX_NACA_WALL_GRADIENT_WLS": "1",
        "CFDX_NACA_CONVERGED_FORCE_CSV": str(output / "force.csv"),
        "CFDX_NACA_PHYSICAL_HISTORY_CSV": str(output / "history.csv"),
        "CFDX_NACA_FACE_PROFILES_CSV": str(output / "profiles.csv"),
        "CFDX_NACA_CHECKPOINT_PAIRED_MANIFEST": str(manifest),
        "CFDX_NACA_CHECKPOINT_EVERY": str(checkpoint_every),
    })
    if stagnation_diagnostic_only:
        env["CFDX_NACA_STAGNATION_DIAGNOSTIC_ONLY"] = "1"
    if initial_paired_checkpoint is not None:
        env["CFDX_NACA_IMPORT_PAIRED_MANIFEST"] = str(
            initial_paired_checkpoint.resolve())
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
            "stagnation_diagnostic_only": stagnation_diagnostic_only,
            "nonlinear_tolerance": 1e-8, "max_iter": maxit,
            "runtime_seconds": round(elapsed, 2),
            "exit_code": proc.returncode,
            "certificate_present": (output / "force.csv").exists(),
            "history_present": (output / "history.csv").exists(),
            "restart_state_used_as_initial_guess": (
                str(initial_state.resolve()) if initial_state else None),
            "initial_paired_manifest": (
                str(initial_paired_checkpoint.resolve())
                if initial_paired_checkpoint else None),
            "last_paired_checkpoint_manifest": (
                str(manifest.resolve()) if manifest.is_file() else None),
            "restart_face_flux_input": (
                str(initial_face_flux.resolve()) if initial_face_flux else None),
            "restart_fidelity": (
                "PAIRED_ACCEPTED_RHIE_CHOW_FACE_FLUX"
                if initial_paired_checkpoint
                else "LEGACY_PAIRED_FACE_FLUX_UNVERIFIED_ITERATION"
                if initial_face_flux
                else "CELL_STATE_ONLY_FACE_FLUX_RECONSTRUCTED"
                if initial_state else "COLD_START"),
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
        "--initial-paired-checkpoint", type=Path,
        help="atomic manifest for same-mesh accepted U,p,Rhie-Chow flux")
    parser.add_argument("--checkpoint-every", type=int, default=300)
    parser.add_argument(
        "--stagnation-diagnostic-only", action="store_true",
        help="do not stop at stagnation; final physical criteria unchanged")
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
        args.initial_face_flux, args.export_converged_face_flux,
        args.initial_paired_checkpoint, args.checkpoint_every,
        args.stagnation_diagnostic_only))
