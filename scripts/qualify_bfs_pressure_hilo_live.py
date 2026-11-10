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

_PREFIX = re.compile(r"^(PRESSURE_HILO_LIVE|PRESSURE_HILO_FACES_LIVE|PRESSURE_HILO_RHIE_CHOW_BRIDGE)\s+(.*)$")
_KV = re.compile(r"([a-zA-Z_][a-zA-Z_0-9]*)=([^\s]+)")


def parse_log(log: str, threshold: float = 1e-11) -> dict:
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("pressure threshold must be finite and positive")
    observations = {"PRESSURE_HILO_LIVE": [], "PRESSURE_HILO_FACES_LIVE": [],
                    "PRESSURE_HILO_RHIE_CHOW_BRIDGE": []}
    for line in log.splitlines():
        match = _PREFIX.match(line)
        if not match:
            continue
        row = dict(_KV.findall(match.group(2)))
        for field in ("iter", "corrector"):
            if field not in row or not row[field].isdigit():
                raise ValueError(f"missing numerical {field} in {match.group(1)}")
        if match.group(1) == "PRESSURE_HILO_LIVE":
            required = ("original_double_relative", "final_split_relative", "final_fp64_relative")
        elif match.group(1) == "PRESSURE_HILO_FACES_LIVE":
            required = ("before_relative", "after_relative", "original_csr_after_relative")
        elif row.get("status") == "UNSUPPORTED_ALGORITHM_OR_SKEW":
            required = ()
        else:
            required = ("rhs_base_flux_relative", "implicit_linearized_continuity_relative",
                        "implicit_operator_relative", "full_rc_compensated_continuity_relative",
                        "full_rc_fp64_continuity_relative", "full_rc_face_rounding_l2_over_rhs")
        for field in required:
            if field not in row or not math.isfinite(float(row[field])) or float(row[field]) < 0:
                raise ValueError(f"nonfinite/missing residual {field}")
        gap_fields = (
            "deferred_divergence_relative", "deferred_face_l2_over_rhs",
            "worst_deferred_face", "worst_deferred_cell", "worst_full_divergence_cell",
            "worst_full_divergence_x", "worst_full_divergence_y",
            "worst_full_divergence_z")
        if any(key in row for key in gap_fields):
            if match.group(1) != "PRESSURE_HILO_RHIE_CHOW_BRIDGE":
                raise ValueError("unexpected spatial diagnostic outside Rhie-Chow")
            for key in gap_fields:
                if key not in row or not math.isfinite(float(row[key])):
                    raise ValueError(f"nonfinite/missing spatial bridge field {key}")
            for key in ("deferred_divergence_relative", "deferred_face_l2_over_rhs"):
                if float(row[key]) < 0:
                    raise ValueError(f"negative bridge residual {key}")
            for key in ("worst_deferred_face", "worst_deferred_cell",
                        "worst_full_divergence_cell"):
                if not row[key].isdigit():
                    raise ValueError(f"invalid bridge index {key}")
        observations[match.group(1)].append(row)
    split = observations["PRESSURE_HILO_LIVE"]
    faces = observations["PRESSURE_HILO_FACES_LIVE"]
    bridge = observations["PRESSURE_HILO_RHIE_CHOW_BRIDGE"]
    if not split or len(split) != len(faces):
        raise ValueError("missing or unpaired live CSR/face diagnostics")
    if bridge and len(bridge) != len(faces):
        raise ValueError("unpaired Rhie-Chow bridge diagnostics")
    comparisons = []
    for i, (x, y) in enumerate(zip(split, faces)):
        if (x["iter"], x["corrector"]) != (y["iter"], y["corrector"]):
            raise ValueError("CSR and face diagnostics refer to different FV states")
        if x.get("physical_fv_qualified") != "0" or y.get("full_rhie_chow_qualified") != "0":
            raise ValueError("unexpected physical qualification in algebraic diagnostic")
        csr_error = float(y["original_csr_after_relative"])
        face_error = float(y["after_relative"])
        extra = {"available": False, "status": "MISSING"}
        if bridge:
            z = bridge[i]
            if (z["iter"], z["corrector"]) != (x["iter"], x["corrector"]):
                raise ValueError("Rhie-Chow bridge refers to a different FV state")
            if z.get("full_rc_qualified") != "0" or z.get("pressure_state_unchanged") != "1":
                raise ValueError("Rhie-Chow bridge cannot qualify or modify physical state")
            if z.get("status") == "UNSUPPORTED_ALGORITHM_OR_SKEW":
                extra = {"available": False, "status": "UNSUPPORTED_ALGORITHM_OR_SKEW"}
            else:
                extra = {"available": True, "status": "DIAGNOSTIC_ONLY",
                         "rhs_base_flux_relative": float(z["rhs_base_flux_relative"]),
                         "implicit_linearized_continuity_relative": float(
                             z["implicit_linearized_continuity_relative"]),
                         "full_rc_compensated_continuity_relative": float(
                             z["full_rc_compensated_continuity_relative"]),
                         "full_rc_fp64_continuity_relative": float(
                             z["full_rc_fp64_continuity_relative"]),
                         "rhs_base_contract_pass": z.get("rhs_base_contract_pass") == "1"}
                if "deferred_divergence_relative" in z:
                    extra["deferred_divergence_relative"] = float(
                        z["deferred_divergence_relative"])
                    extra["deferred_face_l2_over_rhs"] = float(
                        z["deferred_face_l2_over_rhs"])
                    for index in ("worst_deferred_face", "worst_deferred_cell",
                                  "worst_full_divergence_cell"):
                        extra[index] = int(z[index])
                    for axis in "xyz":
                        name = "worst_full_divergence_" + axis
                        extra[name] = float(z[name])
        comparisons.append({
            "iteration": int(x["iter"]), "corrector": int(x["corrector"]),
            "rhie_chow_bridge": extra,
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
        "n_rhie_chow_bridge_diagnostics": len(bridge),
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
