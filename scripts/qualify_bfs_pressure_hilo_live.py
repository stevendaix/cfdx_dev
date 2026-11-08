"""Fail-closed accounting for live production BFS pressure high/low diagnostics.

This is *not* a CFD qualification. It measures whether the rounded CSR and
physical implicit-face operators agree and are both solved, and NEVER promotes
an algebraic certificate to full Rhie-Chow/velocity/SST convergence.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

_PREFIX = re.compile(r"^(PRESSURE_HILO_LIVE|PRESSURE_HILO_FACES_LIVE)\s+(.*)$")
_KV = re.compile(r"([a-zA-Z_][a-zA-Z_0-9]*)=([^\s]+)")


def parse_log(log: str, threshold: float = 1e-11) -> dict:
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("pressure threshold must be finite and positive")
    observations = {"PRESSURE_HILO_LIVE": [], "PRESSURE_HILO_FACES_LIVE": []}
    for line in log.splitlines():
        match = _PREFIX.match(line)
        if not match:
            continue
        row = dict(_KV.findall(match.group(2)))
        for field in ("iter", "corrector"):
            if field not in row or not row[field].isdigit():
                raise ValueError(f"missing numerical {field} in {match.group(1)}")
        for field in (("original_double_relative", "final_split_relative", "final_fp64_relative")
                      if match.group(1) == "PRESSURE_HILO_LIVE"
                      else ("before_relative", "after_relative", "original_csr_after_relative")):
            if field not in row or not math.isfinite(float(row[field])) or float(row[field]) < 0:
                raise ValueError(f"nonfinite/missing residual {field}")
        observations[match.group(1)].append(row)
    split = observations["PRESSURE_HILO_LIVE"]
    faces = observations["PRESSURE_HILO_FACES_LIVE"]
    if not split or len(split) != len(faces):
        raise ValueError("missing or unpaired live CSR/face diagnostics")
    comparisons = []
    for x, y in zip(split, faces):
        if (x["iter"], x["corrector"]) != (y["iter"], y["corrector"]):
            raise ValueError("CSR and face diagnostics refer to different FV states")
        if x.get("physical_fv_qualified") != "0" or y.get("full_rhie_chow_qualified") != "0":
            raise ValueError("unexpected physical qualification in algebraic diagnostic")
        csr_error = float(y["original_csr_after_relative"])
        face_error = float(y["after_relative"])
        comparisons.append({
            "iteration": int(x["iter"]), "corrector": int(x["corrector"]),
            "csr_split_relative": float(x["final_split_relative"]),
            "face_implicit_relative": face_error,
            "csr_after_face_refinement_relative": csr_error,
            "csr_split_pass": float(x["final_split_relative"]) <= threshold,
            "face_implicit_pass": face_error <= threshold,
            "both_face_and_csr_after_face_refinement_pass": max(csr_error, face_error) <= threshold,
            "production_fp64_pass": float(x["original_double_relative"]) <= threshold,
        })
    return {
        "schema": "cfdx-bfs-pressure-hilo-live/v1", "threshold": threshold,
        "n_observations": len(comparisons), "comparisons": comparisons,
        "physical_FV_Rhie_Chow_qualified": False,
        "SST_yplus_Cf_reattachment_qualified": False,
        "status": "DIAGNOSTIC_ONLY_NOT_QUALIFIED",
        "complete_flux_gate": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=1e-11)
    args = parser.parse_args()
    try:
        report = parse_log(args.log.read_text(encoding="utf-8", errors="replace"),
                           args.threshold)
    except (OSError, ValueError) as exc:
        report = {"schema": "cfdx-bfs-pressure-hilo-live/v1", "status": "INCOMPLETE",
                  "error": str(exc), "physical_FV_Rhie_Chow_qualified": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 2  # Diagnostic data can never certify a full CFD simulation.


if __name__ == "__main__":
    raise SystemExit(main())
