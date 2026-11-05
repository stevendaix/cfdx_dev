"""Reproducible A/B for dimensionless SST gates. No physical-qualification claim."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path

import h5py
from generate_th03_openfoam_blockmesh import make_structured_mesh, write_case

PFAIL = re.compile(
    r"pressure-correction solve did not converge .*?relative=([+\d.eE-]+)"
)
SFAIL = re.compile(
    r"SST update rejected:.*?iterations=(\d+).*?k_residual=([+\d.eE-]+), omega_residual=([+\d.eE-]+)"
)


def run(solver: Path, output: Path, timeout: int, iterations: int, inner: int):
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for label, mode, pressure_target in (
        ("legacy_diagnostic", False, 1.0e-9),
        ("relative_diagnostic", True, 1.0e-9),
        ("relative_strict_pressure", True, 1.0e-11),
    ):
        h5_path = output / f"{label}.cfdx.h5"
        write_case(
            h5_path,
            make_structured_mesh("smoke"),
            mesh_source="OpenFOAM six-block reconstruction (smoke); not benchmark-qualified",
            max_iterations=iterations,
            pressure_krylov="cg",
            pressure_preconditioner="native_amg",
            pressure_tolerance=pressure_target,
        )
        with h5py.File(h5_path, "r+") as handle:
            setup = json.loads(handle.attrs["case_setup_json"])
            settings = setup["numerics"]["raw_settings"]
            settings["sst_transport_max_iterations"] = str(inner)
            settings["sst_picard_damping"] = "0.2"
            settings["sst_scale_aware_convergence"] = "true" if mode else "false"
            settings["pressure_face_interpolation"] = "projected"
            handle.attrs["case_setup_json"] = json.dumps(setup)
        command = [
            str(solver.resolve()),
            "--mesh",
            str(h5_path),
            "--output-dir",
            str(output / f"out_{label}"),
            "--iterations",
            str(iterations),
            "--experimental-sst",
            "--benchmark-export",
        ]
        start = time.monotonic()
        timed_out = False
        try:
            completed = subprocess.run(
                command,
                text=True,
                capture_output=True,
                env={**os.environ, "CFDX_SST_DIAGNOSTICS": "1"},
                timeout=timeout,
                check=False,
            )
            log = completed.stdout + completed.stderr
            exit_code = completed.returncode
        except subprocess.TimeoutExpired as error:
            timed_out = True
            log = ((error.stdout or b"") + (error.stderr or b"")).decode(
                "utf-8", "replace"
            )
            exit_code = None
        elapsed = time.monotonic() - start
        (output / f"{label}.log").write_text(log)
        pressure_fail = PFAIL.findall(log)
        sst_fail = SFAIL.findall(log)
        result = {
            "variant": label,
            "pressure_target": pressure_target,
            "dimensionless_sst": mode,
            "sst_inner_budget": inner,
            "return_code": exit_code,
            "timeout": timed_out,
            "elapsed_seconds": round(elapsed, 3),
            "last_pressure_relative_failure": float(pressure_fail[-1])
            if pressure_fail
            else None,
            "last_sst_failure": {
                "inner_iterations": int(sst_fail[-1][0]),
                "k_raw_residual": float(sst_fail[-1][1]),
                "omega_raw_residual": float(sst_fail[-1][2]),
            }
            if sst_fail
            else None,
            "physical_validation": "NOT_QUALIFIED",
            "cf_cp_xr": None,
        }
        results.append(result)
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--iterations", type=int, default=2)
    parser.add_argument("--sst-inner-budget", type=int, default=100)
    a = parser.parse_args()
    print(
        json.dumps(
            run(a.solver, a.output, a.timeout, a.iterations, a.sst_inner_budget),
            indent=2,
        )
    )
