"""Fail-closed NACA A/B campaign: actual solver results, timings, proof paths.

No PASS on a subprocess exit code alone; the physical single-mesh PASS marker,
finite force coefficients and a continuity gate are mandatory.  The final
three-grid GCI remains an independent gate, not inferred from this launcher.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time


_RESULT = re.compile(
    r"NACA0012_RE1000_A0\s+level=\d+\s+Cd=(\S+)\s+Cl=(\S+)"
    r"\s+Cm=(\S+)\s+Cd_pressure=(\S+)\s+Cd_viscous=(\S+)"
    r"\s+continuity=(\S+)\s+iterations=(\d+)"
)
_MODES = {
    "baseline": {},
    "skew": {"CFDX_NACA_SKEW_GRADIENT": "1"},
    "line-search": {"CFDX_NACA_PRESSURE_LINE_SEARCH": "1"},
    "affine-flux": {
        "CFDX_NACA_PRESSURE_LINE_SEARCH": "1",
        "CFDX_NACA_AFFINE_FLUX_TRIALS": "1",
    },
    "tvd": {"CFDX_NACA_CONVECTION_TVD": "1"},
    "wall-face": {"CFDX_NACA_WALL_TRACTION_FACE": "1"},
    "split-boundary": {},
    "split-wall-face": {"CFDX_NACA_WALL_TRACTION_FACE": "1"},
}


def parse_output(output: str, returncode: int) -> dict:
    match = list(_RESULT.finditer(output))
    if returncode != 0 or "NACA0012_SINGLE: PASS" not in output or len(match) != 1:
        raise ValueError("missing verified physical single-mesh success")
    x = match[0]
    values = [float(x.group(i)) for i in range(1, 7)]
    if not all(math.isfinite(v) for v in values):
        raise ValueError("nonfinite force or continuity")
    cd, cl, cm, cdp, cdv, continuity = values
    if not (continuity < 1e-7) or abs(cl) > 0.02 or cdp <= 0:
        raise ValueError("NACA physical gate rejected result")
    if not math.isclose(cd, cdp + cdv, rel_tol=1e-7, abs_tol=1e-9):
        raise ValueError("force decomposition is not conservative")
    return {
        "cd": cd, "cl": cl, "cm": cm, "cd_pressure": cdp,
        "cd_viscous": cdv, "continuity": continuity,
        "iterations": int(x.group(7)),
    }


def run_one(executable: Path, mesh: Path, mode: str, seconds: int, log: Path) -> dict:
    env = os.environ.copy()
    for key in {
        "CFDX_NACA_SKEW_GRADIENT", "CFDX_NACA_PRESSURE_LINE_SEARCH",
        "CFDX_NACA_AFFINE_FLUX_TRIALS", "CFDX_NACA_CONVECTION_TVD",
        "CFDX_NACA_WALL_TRACTION_FACE", "CFDX_NACA_FACE_PROFILES_CSV",
    }:
        env.pop(key, None)
    env.update(_MODES[mode])
    if mode in ("wall-face", "split-wall-face"):
        env["CFDX_NACA_FACE_PROFILES_CSV"] = str(log.with_suffix(".csv"))
    start = time.monotonic()
    output, returncode, reason = "", None, ""
    try:
        completed = subprocess.run(
            [str(executable), "--single", str(mesh)],
            env=env, capture_output=True, text=True, timeout=seconds, check=False,
        )
        output = completed.stdout + "\n" + completed.stderr
        returncode = completed.returncode
    except subprocess.TimeoutExpired as exc:
        data = exc.stdout or b""
        error = exc.stderr or b""
        output = (
            data.decode(errors="replace") if isinstance(data, bytes) else data
        ) + "\n" + (
            error.decode(errors="replace") if isinstance(error, bytes) else error
        )
        reason = "TIMEOUT"
    elapsed = time.monotonic() - start
    log.write_text(output, encoding="utf-8")
    result = {"elapsed_s": elapsed, "log": str(log), "status": "FAIL"}
    if reason:
        result["reason"] = reason
    else:
        try:
            result.update(parse_output(output, returncode))
            result["status"] = "PASS"
        except (ValueError, TypeError) as exc:
            result["reason"] = str(exc)
            result["exit_code"] = returncode
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--mesh-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--modes", nargs="+", choices=tuple(_MODES), default=["baseline"])
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    if args.timeout <= 0 or not args.executable.is_file():
        parser.error("executable must exist and timeout be positive")
    args.output.mkdir(parents=True, exist_ok=True)
    campaign = {"status": "PASS", "modes": {}, "all_three_meshes_required": True}
    for mode in args.modes:
        trials = {}
        for n in (64, 128, 256):
            suffix = "_split" if mode in ("split-boundary", "split-wall-face") else ""
            mesh = args.mesh_dir / f"naca0012_laminar_n{n}{suffix}.h5"
            log = args.output / f"{mode}_n{n}.log"
            if mesh.is_file():
                trials[str(n)] = run_one(args.executable, mesh, mode, args.timeout, log)
            else:
                trials[str(n)] = {"status": "NOT_RUN", "reason": "mesh missing"}
            if trials[str(n)]["status"] != "PASS":
                campaign["status"] = "FAIL"
        campaign["modes"][mode] = trials
        # Exact input format of the existing fail-closed three-grid GCI tool.
        # It must independently reject any incomplete/non-asymptotic family.
        all_pass = all(trials[str(n)]["status"] == "PASS" for n in (64,128,256))
        gci_input = {
            "three_grid_solver_gate": "PASS" if all_pass else "FAIL",
            "results": [
                {
                    "N": n,
                    "status": trials[str(n)]["status"],
                    "metrics": {
                        "Cd": trials[str(n)].get("cd"),
                        "Cl": trials[str(n)].get("cl"),
                        "continuity": trials[str(n)].get("continuity"),
                        "Cd_pressure": trials[str(n)].get("cd_pressure"),
                        "Cd_viscous": trials[str(n)].get("cd_viscous"),
                    },
                } for n in (64, 128, 256)
            ],
        }
        (args.output / f"gci_input_{mode}.json").write_text(
            json.dumps(gci_input, indent=2) + "\n", encoding="utf-8"
        )
    (args.output / "result.json").write_text(
        json.dumps(campaign, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0 if campaign["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
