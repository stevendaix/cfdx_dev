"""Fail-closed, artifact-backed physical qualification on three refined meshes.

Only independently retained solver-result JSON can be promoted. This is NOT a
solver, result extractor, CFD reference database, or replacement for CTest.
See docs/validation/THREE_MESH_PHYSICAL_QUALIFICATION.md for the schema and
exact meaning of normalized residual and conservation measurements.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

SHA256 = re.compile(r"[0-9a-f]{64}\Z")
COMMIT = re.compile(r"[0-9a-f]{40}\Z")
PHYSICS = {
    "laminar": ("mass_global", "mass_local"),
    "naca_laminar": ("mass_global", "mass_local"),
    "naca_sst": ("mass_global", "mass_local"),
    "sst": ("mass_global", "mass_local"),
    "thermal": ("energy_global", "energy_local"),
}
NACA_QOIS = {"cd_pressure", "cd_viscous", "cd_total", "cl"}


class EvidenceError(ValueError):
    """An unqualified or incomplete evidence item cannot be promoted."""


def strict_json(path: Path) -> dict[str, Any]:
    def reject_constant(value: str) -> None:
        raise EvidenceError(f"non-finite JSON value: {value}")

    obj = json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)
    if not isinstance(obj, dict):
        raise EvidenceError(f"{path}: expected a JSON object")
    return obj


def positive(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise EvidenceError(f"{label}: expected positive finite number")
    f = float(value)
    if not math.isfinite(f) or f <= 0.0:
        raise EvidenceError(f"{label}: expected positive finite number")
    return f


def nonnegative(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise EvidenceError(f"{label}: expected non-negative finite number")
    f = float(value)
    if not math.isfinite(f) or f < 0.0:
        raise EvidenceError(f"{label}: expected non-negative finite number")
    return f


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise EvidenceError(f"{label}: expected finite number")
    f = float(value)
    if not math.isfinite(f):
        raise EvidenceError(f"{label}: expected finite number")
    return f


def integer(value: Any, label: str) -> int:
    if type(value) is not int or value <= 0:
        raise EvidenceError(f"{label}: expected strictly positive integer")
    return value


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verified_file(root: Path, obj: Any, label: str) -> tuple[Path, str]:
    if not isinstance(obj, dict):
        raise EvidenceError(f"{label}: path and sha256 are required")
    rel = obj.get("path")
    sha = obj.get("sha256")
    if not isinstance(rel, str) or not rel or not isinstance(sha, str) or not SHA256.fullmatch(sha):
        raise EvidenceError(f"{label}: invalid path/sha256")
    if Path(rel).is_absolute():
        raise EvidenceError(f"{label}: absolute paths are forbidden")
    path = (root / rel).resolve()
    if not path.is_relative_to(root) or not path.is_file() or path.stat().st_size == 0:
        raise EvidenceError(f"{label}: missing, empty or out-of-tree artifact: {rel}")
    if digest(path) != sha:
        raise EvidenceError(f"{label}: SHA-256 mismatch: {rel}")
    return path, sha


def _metric_order(
    values: list[float], ratio_coarse: float, ratio_fine: float, name: str
) -> tuple[float | str, float]:
    # Solve the actual *non-uniform-grid* three-level Richardson equation.
    # Assuming equal ratios and applying log(delta1/delta2)/log(r) to unequal
    # grids yields a biased order (and therefore a biased uncertainty).
    coarse, medium, fine = values
    d1 = coarse - medium
    d2 = medium - fine
    if d1 == 0.0 and d2 == 0.0:
        return "mesh_independent", 0.0
    if d1 * d2 <= 0.0 or abs(d1) <= abs(d2):
        raise EvidenceError(f"{name}: non-monotone or non-convergent mesh series")
    target = abs(d1 / d2)
    limit_as_p_goes_to_zero = math.log(ratio_coarse) / math.log(ratio_fine)
    if target <= limit_as_p_goes_to_zero:
        raise EvidenceError(f"{name}: no positive observed-order solution")

    def ratio_at_order(order: float) -> float:
        return math.expm1(order * math.log(ratio_coarse)) / (
            -math.expm1(-order * math.log(ratio_fine))
        )

    low, high = 1.e-8, 12.0
    if ratio_at_order(high) < target:
        raise EvidenceError(f"{name}: implausibly high observed order")
    for _ in range(80):
        mid = (low + high) / 2.0
        if ratio_at_order(mid) < target:
            low = mid
        else:
            high = mid
    order = (low + high) / 2.0
    correction = math.expm1(order * math.log(ratio_fine))
    if not math.isfinite(order) or correction <= 0.0:
        raise EvidenceError(f"{name}: invalid observed order")
    return order, abs(d2) / correction


def qualify(case: dict[str, Any], root: Path, expected_source_sha: str) -> dict[str, Any]:
    if case.get("schema_version") != "1.0":
        raise EvidenceError("schema_version must be 1.0")
    case_id = case.get("case_id")
    physics = case.get("physics")
    sha = case.get("source_sha")
    if not isinstance(case_id, str) or not case_id.strip():
        raise EvidenceError("case_id is required")
    if physics not in PHYSICS:
        raise EvidenceError("unknown physics family")
    if not isinstance(sha, str) or not COMMIT.fullmatch(sha) or sha != expected_source_sha:
        raise EvidenceError("source SHA missing or does not match expected exact revision")
    limits = case.get("criteria")
    references = case.get("references")
    levels = case.get("levels")
    if not isinstance(limits, dict) or not isinstance(references, dict) or not references:
        raise EvidenceError("numerical criteria and independent references are mandatory")
    if not isinstance(levels, list) or len(levels) != 3:
        raise EvidenceError("exactly three mesh levels are required")
    if physics.startswith("naca") and not NACA_QOIS.issubset(references):
        raise EvidenceError("NACA needs separate pressure, viscous, total drag and lift references")

    max_res = positive(limits.get("max_true_residual"), "max_true_residual")
    balance_fields = PHYSICS[physics]
    max_balance = {
        key: positive(limits.get("max_" + key), "max_" + key)
        for key in balance_fields
    }
    min_ratio = positive(limits.get("min_refinement_ratio"), "min_refinement_ratio")
    min_order = positive(limits.get("min_observed_order"), "min_observed_order")
    max_gci = positive(limits.get("max_gci"), "max_gci")
    max_ratio_difference = nonnegative(
        limits.get("max_refinement_ratio_difference"), "max_refinement_ratio_difference"
    )
    if min_ratio <= 1.0 or max_ratio_difference > 0.10:
        raise EvidenceError("invalid refinement-ratio policy")
    cd_closure = None
    if physics.startswith("naca"):
        cd_closure = positive(limits.get("max_drag_decomposition_error"),
                              "max_drag_decomposition_error")

    reference_values: dict[str, tuple[float, float, float]] = {}
    for name, value in references.items():
        if not isinstance(name, str) or not name or not isinstance(value, dict):
            raise EvidenceError("reference metrics must be named objects")
        reference_values[name] = (
            finite(value.get("value"), f"reference {name}"),
            positive(value.get("scale"), f"reference scale {name}"),
            positive(value.get("max_relative_error"), f"reference tolerance {name}"),
        )

    measured: dict[str, list[float]] = {name: [] for name in reference_values}
    hs: list[float] = []
    cells: list[int] = []
    artifacts: list[dict[str, str]] = []
    for idx, level in enumerate(levels):
        label = f"level[{idx}]"
        if not isinstance(level, dict):
            raise EvidenceError(f"{label}: missing level")
        mesh_path, mesh_sha = verified_file(root, level.get("mesh"), f"{label}.mesh")
        result_path, result_sha = verified_file(root, level.get("result"), f"{label}.result")
        result = strict_json(result_path)
        if result.get("source_sha") != sha or result.get("case_id") != case_id or result.get("physics") != physics:
            raise EvidenceError(f"{label}: wrong source revision or case identity")
        if result.get("mesh_sha256") != mesh_sha:
            raise EvidenceError(f"{label}: result does not reference its input mesh")
        if result.get("converged") is not True or result.get("status") != "COMPLETED":
            raise EvidenceError(f"{label}: solver did not converge and complete")
        if nonnegative(result.get("true_residual"), f"{label}.true_residual") > max_res:
            raise EvidenceError(f"{label}: true linear residual too large")
        balances = result.get("balances")
        if not isinstance(balances, dict):
            raise EvidenceError(f"{label}: missing normalized conservation balances")
        for key in balance_fields:
            if nonnegative(balances.get(key), f"{label}.{key}") > max_balance[key]:
                raise EvidenceError(f"{label}: {key} conservation gate failed")
        n = integer(result.get("cells"), f"{label}.cells")
        h = positive(result.get("h"), f"{label}.h")
        if n < 2:
            raise EvidenceError(f"{label}: mesh contains too few cells")
        observed = result.get("qoi")
        if not isinstance(observed, dict) or set(observed) != set(reference_values):
            raise EvidenceError(f"{label}: missing/unexpected QoI; every reference must be evaluated")
        for name in measured:
            measured[name].append(finite(observed[name], f"{label}.{name}"))
        if cd_closure is not None and abs(
            observed["cd_total"] - observed["cd_pressure"] - observed["cd_viscous"]
        ) > cd_closure:
            raise EvidenceError(f"{label}: pressure/viscous drag does not sum to total")
        hs.append(h)
        cells.append(n)
        artifacts.append({"mesh_path": str(mesh_path.relative_to(root)),
                          "mesh_sha256": mesh_sha,
                          "result_path": str(result_path.relative_to(root)),
                          "result_sha256": result_sha})

    if len({a["mesh_sha256"] for a in artifacts}) != 3 or len({a["result_sha256"] for a in artifacts}) != 3:
        raise EvidenceError("duplicate meshes or reused run evidence")
    if not (cells[0] < cells[1] < cells[2] and hs[0] > hs[1] > hs[2]):
        raise EvidenceError("cell counts and characteristic h must refine strictly")
    r0, r1 = hs[0] / hs[1], hs[1] / hs[2]
    if min(r0, r1) < min_ratio or abs(r0 / r1 - 1.0) > max_ratio_difference:
        raise EvidenceError("non-uniform or insufficient mesh refinement ratio")

    metrics: dict[str, Any] = {}
    for name, values in measured.items():
        reference, scale, max_error = reference_values[name]
        order, absolute_gci = _metric_order(values, r0, r1, name)
        if order != "mesh_independent" and order < min_order:
            raise EvidenceError(f"{name}: order {order:.3g} below {min_order}")
        # GCI = Fs * |f_f - f_m| / (r^p - 1); Fs=1.25 for 3 grids.
        gci = 1.25 * absolute_gci / scale
        if gci > max_gci:
            raise EvidenceError(f"{name}: GCI {gci:.3g} above {max_gci}")
        error = abs(values[2] - reference) / scale
        if error > max_error:
            raise EvidenceError(f"{name}: physical reference error {error:.3g} above {max_error}")
        metrics[name] = {
            "values_coarse_to_fine": values,
            "reference": reference,
            "normalization_scale": scale,
            "relative_reference_error": error,
            "observed_order": order,
            "gci_normalized": gci,
        }

    return {
        "schema_version": "1.0",
        "status": "PASS",
        "case_id": case_id,
        "physics": physics,
        "source_sha": sha,
        "mesh_cells": cells,
        "mesh_h": hs,
        "metrics": metrics,
        "artifacts": artifacts,
        "meaning": "artifact-backed gate; external oracle provenance remains the campaign owner's responsibility",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--expected-source-sha", required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    result: dict[str, Any]
    try:
        if not COMMIT.fullmatch(args.expected_source_sha):
            raise EvidenceError("expected exact full 40-character git SHA")
        path = args.input.resolve()
        result = qualify(strict_json(path), path.parent, args.expected_source_sha)
    except (EvidenceError, ValueError, OSError, json.JSONDecodeError, OverflowError) as exc:
        result = {"status": "FAIL", "reason": str(exc)}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(f"CFDX_THREE_MESH_QUALIFICATION: {result['status']}")
    if result["status"] != "PASS":
        print(result["reason"])
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
