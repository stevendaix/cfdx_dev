"""Fail-closed postprocessor for a persisted CFDX-CONVERGENCE history.

This is an *independent audit* of logged evidence, not an alternative solver
stopping criterion and not a proof of spatial convergence. Thresholds MUST be
supplied by the caller from the original physical qualification contract.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def finite_nonnegative(entry: dict[str, Any], key: str) -> float:
    value = entry.get(key)
    if type(value) not in (float, int) or not math.isfinite(value) or value < 0:
        raise ValueError(f"missing, non-finite or negative {key}")
    return float(value)


def audit_history(
    history: dict[str, Any], thresholds: dict[str, float], min_full_order: int
) -> dict[str, Any]:
    if history.get("format") != "CFDX-CONVERGENCE" or history.get("schema_version") != 1:
        raise ValueError("unsupported or absent CFDX-CONVERGENCE format")
    if type(history.get("converged")) is not bool:
        raise ValueError("missing convergence verdict")
    if type(history.get("iterations")) is not int or history["iterations"] < 1:
        raise ValueError("invalid iteration count")
    entries = history.get("history")
    if not isinstance(entries, list) or not entries:
        raise ValueError("missing nonlinear history")
    if type(min_full_order) is not int or min_full_order < 2:
        raise ValueError("at least two full-order iterations required")
    for value in thresholds.values():
        if not math.isfinite(value) or value <= 0:
            raise ValueError("thresholds must be finite and positive")

    last_iter = 0
    full_order_tail = 0
    for sample in entries:
        if not isinstance(sample, dict):
            raise ValueError("non-object history entry")
        iteration = sample.get("iteration")
        if type(iteration) is not int or iteration <= last_iter:
            raise ValueError("non-monotone or invalid iteration numbers")
        last_iter = iteration
        fraction = finite_nonnegative(sample, "high_order_fraction")
        if fraction > 1:
            raise ValueError("high-order fraction exceeds one")
        full_order = fraction == 1.0
        if full_order:
            full_order_tail += 1
        else:
            full_order_tail = 0
        if type(sample.get("full_order_eligible")) is not bool:
            raise ValueError("missing authoritative full-order eligibility")
        if sample["full_order_eligible"] and not full_order:
            raise ValueError("intermediate order incorrectly marked eligible")
    if last_iter != history["iterations"]:
        raise ValueError("last history iteration != reported iteration")
    final = entries[-1]
    metrics = {
        "momentum_relative": finite_nonnegative(
            final, "momentum_equation_residual_rhs_active_component_max"
        ),
        "pressure_relative": finite_nonnegative(final, "pressure_residual"),
        "continuity_relative": finite_nonnegative(final, "continuity_normalized"),
        "mass_relative": finite_nonnegative(final, "mass_normalized_imbalance"),
    }
    failures: list[str] = []
    if history["converged"] is not True:
        failures.append("solver_not_converged")
    if final["high_order_fraction"] != 1.0:
        failures.append("final_operator_incomplete")
    if not final["full_order_eligible"] or full_order_tail < min_full_order:
        failures.append("full_order_hold_insufficient")
    for name, value in metrics.items():
        if value > thresholds[name]:
            failures.append(f"{name}_above_limit")

    return {
        "qualified": not failures,
        "kind": "nonlinear_iteration_and_conservation_evidence_only",
        "iterations": last_iter,
        "full_order_tail": full_order_tail,
        "metrics": metrics,
        "thresholds": thresholds,
        "failures": failures,
        "spatial_convergence": "not_evaluated",
        "force_accuracy": "not_evaluated",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", type=Path)
    parser.add_argument("--momentum-relative", type=float, required=True)
    parser.add_argument("--pressure-relative", type=float, required=True)
    parser.add_argument("--continuity-relative", type=float, required=True)
    parser.add_argument("--mass-relative", type=float, required=True)
    parser.add_argument("--min-full-order", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    thresholds = {
        "momentum_relative": args.momentum_relative,
        "pressure_relative": args.pressure_relative,
        "continuity_relative": args.continuity_relative,
        "mass_relative": args.mass_relative,
    }
    try:
        with args.history.open(encoding="utf-8") as stream:
            data = json.load(stream, parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
        report = audit_history(data, thresholds, args.min_full_order)
    except (OSError, ValueError, TypeError) as exc:
        report = {"qualified": False, "failures": [f"invalid_evidence: {exc}"]}
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0 if report["qualified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
