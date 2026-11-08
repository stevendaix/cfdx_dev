"""Fail-closed NACA nonlinear acceleration A/B evidence comparison.

The comparator NEVER certifies a CFD solution: a complete official physics
certificate (force/continuity/local-face momentum/spatial) is required in
addition to these time histories. Its purpose is to reject incomparable or
corrupted histories and to report what was actually computed.

Input manifest:
  {"binary":"naca_solver", "mesh":"mesh.h5", "alpha_u":0.6,
   "alpha_p":0.25, "iterations":400, "runs":[
      {"name":"picard", "history":"picard.csv", "method":"OFF"},
      {"name":"anderson", "history":"anderson.csv", "method":"ANDERSON"}
   ]}
Optional per-run expected_hist_sha256, expected_mesh_sha256,
expected_binary_sha256 and interrupted=true (the latter must not certify).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pathlib
import sys
from typing import Any


class EvidenceError(ValueError):
    """A provenance or mathematical gate failed."""


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def positive_number(value: Any, where: str, *, zero: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise EvidenceError(f"{where}: missing or invalid number") from exc
    if not math.isfinite(number) or (number < 0 if zero else number <= 0):
        raise EvidenceError(f"{where}: nonfinite or out-of-range number")
    return number


def strict_float(row: dict[str, str], key: str) -> float:
    if key not in row:
        raise EvidenceError(f"CSV missing required {key}")
    try:
        value = float(row[key])
    except (TypeError, ValueError) as exc:
        raise EvidenceError(f"CSV invalid {key}") from exc
    if not math.isfinite(value):
        raise EvidenceError(f"CSV nonfinite {key}")
    return value


NUMERIC_FIELDS = (
    "continuity_normalized", "momentum_rhs_relative", "mass_normalized_imbalance",
    "mass_local_linf", "face_momentum_x_relative", "face_momentum_y_relative",
    "momentum_x_rhs_relative", "momentum_y_rhs_relative", "momentum_component_max",
    "velocity_change_inf", "pressure_change_inf", "pressure_residual",
    "alpha_u", "alpha_p", "high_order_fraction", "multisecant_relative_step",
    "multisecant_estimated_objective", "nonfinite_mass_faces", "nonfinite_velocities",
    "multisecant_proposed", "multisecant_accepted",
)


def parse_run(root: pathlib.Path, run: dict[str, Any],
              global_config: dict[str, Any], shared_hashes: dict[str, str]) -> dict[str, Any]:
    name = str(run.get("name", ""))
    method = str(run.get("method", ""))
    if not name or method not in {"OFF", "ANDERSON", "NGMRES"}:
        raise EvidenceError("run requires a unique name and method OFF/ANDERSON/NGMRES")
    for k in ("mesh", "binary"):
        expected = run.get(f"expected_{k}_sha256", global_config.get(f"expected_{k}_sha256"))
        if expected is not None and expected != shared_hashes[k]:
            raise EvidenceError(f"{name}: {k} sha256 mismatch")
    history_path = (root / str(run.get("history", ""))).resolve()
    if not history_path.is_file():
        raise EvidenceError(f"{name}: history missing")
    history_hash = sha256(history_path)
    if "expected_hist_sha256" in run and history_hash != run["expected_hist_sha256"]:
        raise EvidenceError(f"{name}: history SHA-256 mismatch")
    with history_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise EvidenceError(f"{name}: empty CSV")
    for i, row in enumerate(rows, 1):
        if strict_float(row, "iteration") != i:
            raise EvidenceError(f"{name}: iteration missing, duplicate or shuffled at {i}")
        for k in NUMERIC_FIELDS:
            strict_float(row, k)
        if row.get("convection_scheme") != "TVD" or \
           abs(strict_float(row, "high_order_fraction")-1.0) > 1e-12:
            raise EvidenceError(f"{name}: final second-order TVD operator not used in row {i}")
        if row.get("algorithm") != "SIMPLE":
            raise EvidenceError(f"{name}: non-SIMPLE algorithm")
        if row.get("multisecant_method") != method:
            raise EvidenceError(f"{name}: method in CSV does not match manifest")
        if abs(strict_float(row, "alpha_u")-global_config["alpha_u"]) > 1e-13 or \
           abs(strict_float(row, "alpha_p")-global_config["alpha_p"]) > 1e-13:
            raise EvidenceError(f"{name}: relaxation changed during comparison")
        if strict_float(row, "nonfinite_mass_faces") != 0 or \
           strict_float(row, "nonfinite_velocities") != 0:
            raise EvidenceError(f"{name}: nonfinite flow state")
        if strict_float(row, "multisecant_accepted") not in {0, 1} or \
           strict_float(row, "multisecant_proposed") not in {0, 1}:
            raise EvidenceError(f"{name}: malformed acceleration flags")
        if method == "OFF" and strict_float(row, "multisecant_accepted") != 0:
            raise EvidenceError(f"{name}: OFF run falsely reports acceleration")
    final = rows[-1]
    total = int(global_config["iterations"])
    complete = len(rows) == total
    if len(rows) > total:
        raise EvidenceError(f"{name}: history longer than requested budget")
    if not complete and not run.get("interrupted", False):
        raise EvidenceError(f"{name}: incomplete history not marked interrupted")
    if complete and run.get("interrupted", False):
        raise EvidenceError(f"{name}: complete history incorrectly marked interrupted")
    extracted = {
        "momentum_global": strict_float(final, "momentum_rhs_relative"),
        "momentum_x": strict_float(final, "momentum_x_rhs_relative"),
        "momentum_y": strict_float(final, "momentum_y_rhs_relative"),
        "face_momentum_x": strict_float(final, "face_momentum_x_relative"),
        "face_momentum_y": strict_float(final, "face_momentum_y_relative"),
        "continuity": strict_float(final, "continuity_normalized"),
        "local_mass": strict_float(final, "mass_local_linf"),
        "mass_imbalance": strict_float(final, "mass_normalized_imbalance"),
        "velocity_change": strict_float(final, "velocity_change_inf"),
    }
    tol = global_config["tolerance"]
    physics_gates = {
        "momentum": max(extracted["momentum_x"], extracted["momentum_y"],
                        extracted["momentum_global"]) <= tol,
        "independent_face_momentum": max(extracted["face_momentum_x"],
                                         extracted["face_momentum_y"]) <= tol,
        "continuity": extracted["continuity"] <= tol,
        "mass_imbalance": extracted["mass_imbalance"] <= tol,
        "local_mass": extracted["local_mass"] <= tol,
        "velocity_change": extracted["velocity_change"] <= tol,
    }
    accepted = sum(int(float(row["multisecant_accepted"])) for row in rows)
    return {
        "name": name, "method": method, "history": str(history_path),
        "sha256": history_hash, "iterations": len(rows), "budget": total,
        "complete": complete, "accelerated_accepted": accepted,
        "last_metrics": extracted, "physics_gates": physics_gates,
        # Deliberate fail-closed: a history has no authoritative certification.
        "qualification": "NOT_CERTIFIED", "requires_official_certificate": True,
    }


def compare(manifest_path: pathlib.Path) -> dict[str, Any]:
    manifest_path = manifest_path.resolve()
    root = manifest_path.parent
    config = json.loads(manifest_path.read_text(encoding="utf-8"))
    config["alpha_u"] = positive_number(config.get("alpha_u"), "alpha_u")
    config["alpha_p"] = positive_number(config.get("alpha_p"), "alpha_p")
    config["tolerance"] = positive_number(config.get("tolerance", 1e-8), "tolerance")
    if config["tolerance"] > 1e-8:
        raise EvidenceError("strict physics tolerance must not be relaxed above 1e-8")
    iterations = positive_number(config.get("iterations"), "iterations")
    if iterations != int(iterations):
        raise EvidenceError("iterations must be integer")
    config["iterations"] = int(iterations)
    shared_hashes = {}
    for name in ("binary", "mesh"):
        path = (root / str(config.get(name, ""))).resolve()
        if not path.is_file():
            raise EvidenceError(f"common {name} missing")
        shared_hashes[name] = sha256(path)
    runs = config.get("runs")
    if not isinstance(runs, list) or not runs:
        raise EvidenceError("manifest needs runs")
    output = [parse_run(root, run, config, shared_hashes) for run in runs]
    if len(set(result["name"] for result in output)) != len(output):
        raise EvidenceError("duplicate run names")
    if not any(result["method"] == "OFF" for result in output):
        raise EvidenceError("baseline Picard OFF missing")
    baseline = next(result for result in output if result["method"] == "OFF")
    for result in output:
        result["relative_to_picard"] = {
            key: (baseline["last_metrics"][key] / value if value != 0 else None)
            for key, value in result["last_metrics"].items()
        }
    return {
        "status": "EVIDENCE_COMPARABLE_NO_PHYSICAL_QUALIFICATION",
        "same_binary_sha256": shared_hashes["binary"],
        "same_mesh_sha256": shared_hashes["mesh"],
        "strict_tolerance": config["tolerance"],
        "controls": {"alpha_u": config["alpha_u"],
                     "alpha_p": config["alpha_p"], "iteration_budget": config["iterations"]},
        "runs": output,
        "warning": "A physical run is never certified using CSV time histories alone.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path)
    args = parser.parse_args()
    try:
        outcome = compare(args.manifest)
    except (EvidenceError, OSError, ValueError, TypeError, KeyError) as exc:
        print(f"MULTISECANT_EVIDENCE_FAIL: {exc}", file=sys.stderr)
        return 2
    target = json.dumps(outcome, indent=2, ensure_ascii=False, sort_keys=True)
    if args.output:
        args.output.write_text(target + "\n", encoding="utf-8")
    print("MULTISECANT_EVIDENCE_VALID; NO PHYSICAL QUALIFICATION")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
