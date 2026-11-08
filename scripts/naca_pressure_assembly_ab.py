"""Production NACA SIMPLE/MUSCL pressure-assembly A/B; fail-closed.

This is a one-grid regression comparison, not a N128/N256/GCI certificate.
The same production executable and mesh are used, varying only the explicitly
opted-in compensated pressure assembly. Physical gates remain in the driver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path

COLUMNS = ("cd", "cdp", "cdv", "cl", "cm")
BASE_ENV = {
    "CFDX_NACA_ALGORITHM": "SIMPLE",
    "CFDX_NACA_CONVECTION_CELL_LIMITED": "1",
    "CFDX_NACA_FARFIELD_PRESSURE_NEUMANN": "1",
    "CFDX_NACA_REQUIRE_KINEMATIC_FACE_FLUX": "1",
    "CFDX_NACA_LINEAR_TOLERANCE": "1e-12",
    "CFDX_NACA_WALL_TRACTION_FACE": "1",
    "CFDX_NACA_WALL_GRADIENT_WLS": "1",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for buf in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(buf)
    return h.hexdigest()


def parse_results(run_dir: Path) -> dict:
    log = (run_dir / "run.log").read_text(encoding="utf-8")
    wall = re.findall(
        r"^NACA_SOLVE_WALL .*accepted_iterations=(\d+) converged=(\d+) algorithm=(\w+)",
        log,
        flags=re.M,
    )
    if len(wall) != 1 or wall[0][1:] != ("1", "SIMPLE"):
        raise ValueError("NACA A/B: solver did not physically converge with SIMPLE")
    if not re.search(
        r"^NACA_FARFIELD_FLUX_AUDIT .*kinematic_match=1 .*physically_converged=1",
        log,
        re.M,
    ):
        raise ValueError("NACA A/B: missing physical farfield flux certificate")
    # The five-argument --single continuation branch reports its own
    # converged force certificate, rather than the legacy two-argument PASS.
    # Require the actual production certificate and its matching iteration.
    if not re.search(
        r"^NACA_CONVERGED_FORCE_CERTIFICATE .*qualification=ONE_GRID_ONLY$", log, re.M
    ):
        raise ValueError("NACA A/B: missing one-grid converged force certificate")
    reported_iterations = re.findall(
        r"^NACA0012_RE1000_A0 .*iterations=(\d+)$", log, re.M
    )
    if reported_iterations != [wall[0][0]]:
        raise ValueError("NACA A/B: NACA result line missing or iteration mismatch")
    with (run_dir / "force.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    if len(rows) != 1 or rows[0].get("converged") != "yes":
        raise ValueError("NACA A/B: missing unique converged force certificate")
    with (run_dir / "history.csv").open(newline="", encoding="utf-8") as stream:
        hist = list(csv.DictReader(stream))
    if (
        len(hist) != int(wall[0][0])
        or not hist
        or any(
            row.get("algorithm") != "SIMPLE"
            or row.get("convection_scheme") != "MUSCL_CELL_LIMITED"
            for row in hist
        )
    ):
        raise ValueError("NACA A/B: incomplete history or scheme mismatch")
    if int(hist[-1]["iteration"]) != len(hist):
        raise ValueError("NACA A/B: final iteration mismatch")
    force = {key: float(rows[0][key]) for key in COLUMNS}
    force["res_mom"] = float(rows[0]["res_mom"])
    force["res_mass"] = float(rows[0]["res_mass"])
    if any(not math.isfinite(value) for value in force.values()):
        raise ValueError("NACA A/B: nonfinite forces or physical residuals")
    if abs(force["cd"] - force["cdp"] - force["cdv"]) > 1e-9:
        raise ValueError("NACA A/B: pressure/viscous drag not closed")
    return {
        "case": rows[0]["case"],
        "h": float(rows[0]["h"]),
        "iterations": len(hist),
        "force": force,
        "final_high_order_fraction": float(hist[-1]["high_order_fraction"]),
        "final_momentum_component_max": float(hist[-1]["momentum_component_max"]),
        "final_face_momentum_y_relative": float(hist[-1]["face_momentum_y_relative"]),
    }


def compare(reference: dict, candidate: dict, absolute_force_limit: float) -> dict:
    if reference["case"] != candidate["case"] or reference["h"] != candidate["h"]:
        raise ValueError("NACA A/B: different mesh refinement/case")
    if (
        reference["final_high_order_fraction"] != 1.0
        or candidate["final_high_order_fraction"] != 1.0
    ):
        raise ValueError("NACA A/B: final order-two convection not active")
    differences = {
        key: candidate["force"][key] - reference["force"][key] for key in COLUMNS
    }
    if not math.isfinite(absolute_force_limit) or absolute_force_limit <= 0:
        raise ValueError("NACA A/B: invalid comparison tolerance")
    return {
        "schema": "cfdx-naca-production-pressure-assembly-ab/v1",
        "status": "ONE_GRID_REGRESSION_PASS"
        if all(abs(v) <= absolute_force_limit for v in differences.values())
        else "ONE_GRID_REGRESSION_FAIL",
        "reference": reference,
        "candidate": candidate,
        "differences": differences,
        "absolute_force_comparison_limit": absolute_force_limit,
        "complete_three_grid_gci": False,
        "same_nonlinear_trajectory": reference["iterations"] == candidate["iterations"],
    }


def one_run(
    executable: Path,
    mesh: Path,
    dest: Path,
    max_iterations: int,
    timeout_seconds: int,
    experimental: bool,
) -> int:
    dest.mkdir(parents=True, exist_ok=True)
    for f in ("run.log", "force.csv", "history.csv", "raw_wall_faces.csv"):
        (dest / f).unlink(missing_ok=True)
    env = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith("CFDX_NACA_") and not k.startswith("CFDX_PRESSURE_")
    }
    env.update(BASE_ENV)
    env["CFDX_NACA_CONVERGED_FORCE_CSV"] = str(dest / "force.csv")
    env["CFDX_NACA_PHYSICAL_HISTORY_CSV"] = str(dest / "history.csv")
    env["CFDX_NACA_RAW_WALL_FACES_CSV"] = str(dest / "raw_wall_faces.csv")
    if experimental:
        env["CFDX_NACA_COMPENSATED_PRESSURE_ASSEMBLY"] = "1"
    with (dest / "run.log").open("w", encoding="utf-8") as stream:
        try:
            proc = subprocess.run(
                [
                    str(executable),
                    "--single",
                    str(mesh),
                    "0.6",
                    "0.25",
                    str(max_iterations),
                ],
                stdout=stream,
                stderr=subprocess.STDOUT,
                env=env,
                timeout=timeout_seconds,
                check=False,
            )
        except subprocess.TimeoutExpired:
            stream.write("\nNACA_AB_TIMEOUT_PHYSICAL_NOT_QUALIFIED\n")
            return 124
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-iterations", default=450, type=int)
    parser.add_argument("--timeout-seconds", default=180, type=int)
    parser.add_argument("--force-comparison-absolute", default=5e-7, type=float)
    parser.add_argument(
        "--resume-verified-baseline",
        action="store_true",
        help="reuse previous same mesh/binary baseline after verifying stored hashes and physical outputs",
    )
    args = parser.parse_args()
    if (
        not args.mesh.is_file()
        or not args.executable.is_file()
        or args.max_iterations <= 0
        or args.timeout_seconds <= 0
    ):
        parser.error("missing executable/mesh or nonpositive budget")
    args.output.mkdir(parents=True, exist_ok=True)
    output = {
        "schema": "cfdx-naca-production-pressure-assembly-ab/v1",
        "mesh_sha256": digest(args.mesh),
        "binary_sha256": digest(args.executable),
        "status": "NOT_QUALIFIED",
        "runs": {},
    }
    old = None
    if args.resume_verified_baseline:
        try:
            old = json.loads(
                (args.output / "comparison.json").read_text(encoding="utf-8")
            )
        except (OSError, ValueError) as exc:
            parser.error(f"missing prior run manifest: {exc}")
        if (
            old.get("mesh_sha256") != output["mesh_sha256"]
            or old.get("binary_sha256") != output["binary_sha256"]
            or old.get("runs", {}).get("baseline", {}).get("return_code") != 0
        ):
            parser.error(
                "baseline reuse forbidden: mesh/binary SHA or completion mismatch"
            )
    for name, experimental in (("baseline", False), ("extended_assembly", True)):
        dest = args.output / name
        if name == "baseline" and old is not None:
            rc = 0  # parse_results below revalidates physical gates on the saved evidence
        else:
            rc = one_run(
                args.executable,
                args.mesh,
                dest,
                args.max_iterations,
                args.timeout_seconds,
                experimental,
            )
        output["runs"][name] = {
            "return_code": rc,
            "pressure_compensated_assembly": experimental,
        }
        if rc != 0:
            output["status"] = "INCOMPLETE_OR_FAILED_RUN"
            break
        try:
            output["runs"][name]["metrics"] = parse_results(dest)
        except (OSError, KeyError, ValueError) as exc:
            output["status"] = "INCOMPLETE_OR_FAILED_RUN"
            output["runs"][name]["error"] = str(exc)
            break
    else:
        try:
            comparison = compare(
                output["runs"]["baseline"]["metrics"],
                output["runs"]["extended_assembly"]["metrics"],
                args.force_comparison_absolute,
            )
            output["comparison"] = comparison
            output["status"] = comparison["status"]
        except ValueError as exc:
            output["status"] = "INVALID_COMPARISON"
            output["comparison_error"] = str(exc)
    (args.output / "comparison.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": output["status"],
                "comparison": output.get("comparison", {}).get("differences", {}),
            },
            indent=2,
        )
    )
    return 0 if output["status"] == "ONE_GRID_REGRESSION_PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
