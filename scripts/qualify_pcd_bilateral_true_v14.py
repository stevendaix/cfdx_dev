#!/usr/bin/env python3
"""Fail-closed consistency audit of real V14 nested PCD/Fp applications.

An accepted projected PRECONDITIONER action does NOT certify outer COUPLED,
SST, or Driver-Seegmiller physics. The original incompatible RHS is retained.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path


class PcdProofError(ValueError):
    pass


def certify(log: str, *, minimum_actions: int = 10, tolerance: float = 1e-8) -> dict:
    if minimum_actions < 1 or not math.isfinite(tolerance) or tolerance <= 0:
        raise PcdProofError("invalid proof controls")
    runs = []
    for line in log.splitlines():
        if not line.startswith("COUPLED_NESTED_PCD_FP "):
            continue
        f = dict(re.findall(r"(\w+)=([^\s]+)", line))
        required = {"bilateral_deflation", "auxiliary_operator", "preconditioner_anchored",
                    "auxiliary_true_relative_rhs", "projected_original_true_relative",
                    "left_incompatible_floor", "original_true_relative_rhs", "accepted",
                    "status", "iterations"}
        if not required <= f.keys():
            raise PcdProofError("incomplete Fp action certificate")
        if f["bilateral_deflation"] != "1" or f["auxiliary_operator"] != "deflated" or f["preconditioner_anchored"] != "1":
            raise PcdProofError("the actual Fp operator/preconditioner does not match bilateral deflation")
        try:
            aux, projected, floor, original = (
                float(f[field]) for field in ("auxiliary_true_relative_rhs",
                    "projected_original_true_relative", "left_incompatible_floor",
                    "original_true_relative_rhs"))
            status = int(f["status"])
            accepted = int(f["accepted"])
            iterations = int(f["iterations"])
        except (ValueError, OverflowError) as exc:
            raise PcdProofError("invalid numeric Fp evidence") from exc
        if (not all(math.isfinite(v) and v >= 0 for v in (aux, projected, floor, original))
            or status != 0 or accepted != 1 or iterations < 0
            or aux > tolerance or projected > tolerance):
            raise PcdProofError("nonconverged or nonfinite bilateral Fp action")
        # Orthogonality: the discarded LEFT-null component remains measurable
        # on the physical/original Fp. Roundoff-safe since logs print 6 digits.
        predicted = math.hypot(floor, projected * math.sqrt(max(0., 1. - floor * floor)))
        if floor > 1 or abs(original - predicted) > 2e-5 * max(1., original):
            raise PcdProofError("true ORIGINAL residual hides incompatible component")
        runs.append({"iterations": iterations, "auxiliary_true": aux,
                     "projected_true": projected, "left_incompatible_floor": floor,
                     "original_true": original})
    if len(runs) < minimum_actions:
        raise PcdProofError(f"only {len(runs)} accepted Fp applications; require {minimum_actions}")
    return {"schema": "CFDX-PCD-BILATERAL-TRUE-V14", "actions": len(runs),
            "max_auxiliary_true": max(r["auxiliary_true"] for r in runs),
            "max_projected_true": max(r["projected_true"] for r in runs),
            "max_original_true": max(r["original_true"] for r in runs),
            "max_incompatible_floor": max(r["left_incompatible_floor"] for r in runs),
            "projected_preconditioner_qualified": True,
            "outer_coupled_converged": False,
            "driver_seegmiller_qualified": False}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-actions", type=int, default=10)
    ns = parser.parse_args()
    try:
        output = certify(ns.log.read_text(errors="replace"), minimum_actions=ns.minimum_actions)
        status = 0
    except (OSError, PcdProofError) as exc:
        output = {"schema": "CFDX-PCD-BILATERAL-TRUE-V14", "projected_preconditioner_qualified": False,
                  "outer_coupled_converged": False, "driver_seegmiller_qualified": False,
                  "reason": str(exc)}
        status = 3
    ns.output.parent.mkdir(parents=True, exist_ok=True)
    ns.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
