"""Fail-closed NACA three-grid diagnostic, not a reference validation.

The solver's --single mode must satisfy its own physical gates. This runner
additionally requires mesh-refinement stabilization; it never converts a
non-converged run, timeout, or missing QoI into PASS.
"""
import argparse
import json
import math
import re
import subprocess
from pathlib import Path


def parse_result(output: str) -> dict:
    if "NACA0012_SINGLE: PASS" not in output:
        raise ValueError("single-mesh solver did not report a validated PASS")
    matches = [
        line
        for line in output.splitlines()
        if line.startswith("NACA0012_RE1000_A0 level=")
    ]
    if len(matches) != 1:
        raise ValueError("expected exactly one NACA coefficient record")
    line = matches[0]
    parsed = {}
    for key in ("Cd", "Cl", "Cd_pressure", "continuity", "iterations"):
        match = re.search(r"(?:^|\s)" + re.escape(key) + r"=([^\s]+)", line)
        if match is None:
            raise ValueError(f"missing {key}")
        value = float(match.group(1))
        if not math.isfinite(value):
            raise ValueError(f"nonfinite {key}")
        parsed[key] = value
    if not (
        parsed["Cd_pressure"] > 0
        and parsed["continuity"] < 1e-7
        and abs(parsed["Cl"]) <= 2e-2
        and parsed["iterations"] >= 1
    ):
        raise ValueError("physical gate failed")
    return parsed


def run_campaign(
    executable: Path, meshes: list[Path], timeout: float, output_dir: Path
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for n, mesh in zip((64, 128, 256), meshes):
        record = {"N": n, "mesh": str(mesh), "status": "FAIL"}
        results.append(record)
        if not mesh.is_file():
            record["reason"] = "missing mesh"
            break
        try:
            proc = subprocess.run(
                [str(executable), "--single", str(mesh)],
                capture_output=True,
                text=True,
                timeout=timeout,
                errors="replace",
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            (output_dir / f"naca_n{n}.stdout.log").write_bytes(exc.stdout or b"")
            (output_dir / f"naca_n{n}.stderr.log").write_bytes(exc.stderr or b"")
            record["reason"] = "timeout"
            break
        (output_dir / f"naca_n{n}.stdout.log").write_text(proc.stdout)
        (output_dir / f"naca_n{n}.stderr.log").write_text(proc.stderr)
        record["exit_code"] = proc.returncode
        if proc.returncode != 0:
            record["reason"] = "nonzero exit"
            break
        try:
            record["metrics"] = parse_result(proc.stdout)
        except ValueError as exc:
            record["reason"] = str(exc)
            break
        record["status"] = "PASS"
    passed = len(results) == 3 and all(r["status"] == "PASS" for r in results)
    if passed:
        d01 = abs(results[1]["metrics"]["Cd"] - results[0]["metrics"]["Cd"])
        d12 = abs(results[2]["metrics"]["Cd"] - results[1]["metrics"]["Cd"])
        passed = d12 <= d01
        if not passed:
            results[-1]["status"] = "FAIL"
            results[-1]["reason"] = "drag not stabilizing under refinement"
    return {
        "three_grid_solver_gate": "PASS" if passed else "FAIL",
        "reference_validation": "NOT_PERFORMED",
        "wall_shear_validation": "NOT_PERFORMED",
        "results": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("mesh_n64", type=Path)
    parser.add_argument("mesh_n128", type=Path)
    parser.add_argument("mesh_n256", type=Path)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--output-dir", type=Path, default=Path("naca_campaign_logs"))
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be positive and finite")
    result = run_campaign(
        args.executable,
        [args.mesh_n64, args.mesh_n128, args.mesh_n256],
        args.timeout,
        args.output_dir,
    )
    (args.output_dir / "qualification_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))
    return 0 if result["three_grid_solver_gate"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
