"""Strict BFS U-p-k-omega proof gate; no physical Driver-Seegmiller PASS by default.

Three levels remain separate:
  1. execution (process and consistent execution.json)
  2. nonlinear (last returned-state V8 evidence and the *resolved* solver plan)
  3. physical comparison (never inferred from this debug grid)

Do not use the case's requested Schur as evidence for its actual execution.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def _reject_constant(value: str):
    raise ValueError(f"non-finite JSON token: {value}")


def _strict_object(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"duplicate JSON field {key!r}")
        result[key] = value
    return result


def strict_load(path: Path) -> dict:
    doc = json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=_reject_constant,
        object_pairs_hook=_strict_object,
    )
    if not isinstance(doc, dict):
        raise TypeError("expected JSON object")
    return doc


def _positive_finite_number(value) -> bool:

    return type(value) in (float, int) and math.isfinite(value) and value >= 0.0


def validate_bfs_evidence(
    output_dir: Path,
    exit_code: int | None,
    *,
    algorithm: str,
    schur: str | None,
    preconditioner: str | None,
    case_sha256: str | None = None,
    case_path: Path | None = None,
) -> dict:
    """Check persisted evidence from this process; never synthesize a PASS.

    A converged debug case can pass the nonlinear *software gate* but cannot
    pass experimental validation. Experimental x_r, Cf and velocity/Reynolds
    stresses require independent samples and an actual reference grid.
    """
    checks = []
    defects = []
    proof = {
        "format": "CFDX-BFS-SCHUR-EVIDENCE/v1",
        "requested_algorithm": algorithm,
        "requested_schur": schur,
        "requested_preconditioner": preconditioner,
        "case_sha256": case_sha256,
        "execution_proven": False,
        "nonlinear_v8_proven": False,
        "physical_validation": "NOT_EVALUATED_NO_EXPERIMENTAL_COMPARISON",
        "measured_reattachment_x_over_h": None,
        "cf_comparison": "NOT_EVALUATED",
        "mean_velocity_comparison": "NOT_EVALUATED",
        "reynolds_stress_comparison": "NOT_EVALUATED",
    }
    if case_path is not None:
        try:
            import h5py

            case_path = Path(case_path)

            actual_sha = hashlib.sha256(case_path.read_bytes()).hexdigest()
            if actual_sha != case_sha256:
                defects.append("case_sha256_mismatch")
            with h5py.File(case_path, "r") as h5:
                setup = json.loads(h5.attrs["case_setup_json"])
                entries = setup["numerics"]["selection"]["entries"]
                families = {x["family"]: x["configuration_key"] for x in entries}
                if len(families) != len(entries):
                    defects.append("duplicate_case_numerical_family")
                if (
                    families.get("pressure_velocity")
                    != f"pressure_velocity.{algorithm}"
                ):
                    defects.append("case_algorithm_mismatch")
                if algorithm == "coupled" and families.get("schur") != f"schur.{schur}":
                    defects.append("case_schur_mismatch")
                if algorithm == "coupled":
                    configured_pc = families.get("preconditioner")
                    expected_pc = (
                        None
                        if preconditioner == "auto"
                        else f"preconditioner.{preconditioner}"
                    )
                    if configured_pc != expected_pc:
                        defects.append("case_preconditioner_mismatch")
                if families.get("convection") not in (
                    "numerics.convection.second_order_upwind",
                    "numerics.convection.tvd.vanleer",
                    "numerics.convection.tvd.vanalbada",
                    "numerics.convection.tvd.mc",
                    "numerics.convection.tvd.superbee",
                ):
                    defects.append("case_second_order_convection_not_selected")
                raw = setup["numerics"]["raw_settings"]
                if raw.get("v8_global_nonlinear_residual") != "true":
                    defects.append("case_v8_certificate_not_requested")
                if raw.get("v7_post_sst_momentum_gate") != "true":
                    defects.append("case_post_sst_momentum_gate_not_requested")
                proof["case_mesh_identity"] = {
                    "topology_hash": str(h5.attrs["topology_hash"]),
                    "geometry_hash": str(h5.attrs["geometry_hash"]),
                    "n_cells": int(h5.attrs["n_cells"]),
                    "n_faces": int(h5.attrs["n_faces"]),
                }
            if not defects:
                checks.append("case_sha256_mesh_and_numerics_verified")
        except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
            defects.append(f"case_provenance_invalid:{exc}")
    if exit_code is None:
        return {
            **proof,
            "checks": checks,
            "defects": defects + ["process_not_completed"],
        }

    execution_path = output_dir / "execution.json"
    convergence_path = output_dir / "convergence.json"
    try:
        summary = strict_load(execution_path)
        history = strict_load(convergence_path)
    except (OSError, ValueError) as exc:
        return {**proof, "checks": [], "defects": [f"evidence_file_invalid:{exc}"]}

    if (
        summary.get("format") != "CFDX-EXECUTION"
        or type(summary.get("process_exit_code")) is not int
        or summary.get("process_exit_code") != exit_code
        or exit_code != 0
    ):
        defects.append("execution_process_not_successful_or_mismatched")
    else:
        checks.append("execution_process_consistent")
        proof["execution_proven"] = True
    if summary.get("converged") is not True:
        defects.append("execution_solver_not_converged")

    if (
        history.get("format") != "CFDX-CONVERGENCE"
        or history.get("converged") is not True
        or type(history.get("iterations")) is not int
        or history["iterations"] < 1
    ):
        defects.append("convergence_history_not_converged")
    else:
        checks.append("convergence_history_converged")

    entries = history.get("history")
    if (
        not isinstance(entries, list)
        or not entries
        or not isinstance(entries[-1], dict)
    ):
        defects.append("last_returned_state_missing")
        return {**proof, "checks": checks, "defects": defects}
    last = entries[-1]
    proof["last_recorded_iteration"] = last.get("iteration")
    if last.get("iteration") != history.get("iterations"):
        defects.append("last_state_not_final_iteration")
    if (
        last.get("full_order_eligible") is not True
        or last.get("high_order_fraction") != 1
    ):
        defects.append("full_second_order_solution_not_evidenced")

    plan = history.get("linear_plan")
    if not isinstance(plan, dict):
        defects.append("actual_linear_plan_missing")
        plan = {}
    if algorithm == "coupled":
        coupled = plan.get("coupled")
        if not isinstance(coupled, dict) or coupled.get("resolved") is not True:
            defects.append("coupled_linear_plan_not_resolved")
            coupled = {}
        else:
            checks.append("coupled_linear_plan_resolved")
        proof["actual_krylov"] = coupled.get("krylov")
        proof["actual_preconditioner"] = coupled.get("preconditioner")
        proof["actual_schur"] = plan.get("resolved_schur")
        if (
            plan.get("schur_resolved") is not True
            or plan.get("requested_schur") != schur
            or plan.get("resolved_schur") != schur
        ):
            defects.append("requested_resolved_schur_mismatch")
        if coupled.get("krylov") not in ("gmres", "fgmres"):
            defects.append("coupled_krylov_missing_or_invalid")
        valid_pcs = ("pcd", "coupled_block_schur", "mgr", "braess_sarazin")
        if coupled.get("preconditioner") not in valid_pcs:
            defects.append("coupled_preconditioner_missing_or_invalid")
        elif preconditioner != "auto" and coupled["preconditioner"] != preconditioner:
            defects.append("requested_resolved_preconditioner_mismatch")
    else:
        if schur is not None or preconditioner is not None:
            defects.append("segregated_run_has_coupled_request")
        if (
            isinstance(plan.get("coupled"), dict)
            and plan["coupled"].get("resolved") is True
        ):
            defects.append("unexpected_coupled_plan_in_segregated_run")
        proof["actual_preconditioner"] = None
        proof["actual_schur"] = None

    v8 = last.get("global_nonlinear_v8")
    if not isinstance(v8, dict):
        defects.append("independent_v8_certificate_missing")
        v8 = {}
    proof["v8"] = v8
    if v8.get("evaluated") is not True:
        defects.append("independent_v8_not_evaluated")
    if v8.get("algorithm") != algorithm:
        defects.append("v8_algorithm_mismatch")
    if algorithm == "coupled" and v8.get("schur") != schur:
        defects.append("v8_schur_mismatch")
    for required in (
        "sst",
        "sst_transport_committed",
        "bounded",
        "all_finite",
        "converged",
    ):
        if v8.get(required) is not True:
            defects.append(f"v8_{required}_not_true")
    for key in (
        "momentum_global_relative",
        "momentum_max_active_component_relative",
        "continuity",
        "k",
        "omega",
        "eddy_viscosity_closure",
        "merit",
    ):
        if not _positive_finite_number(v8.get(key)):
            defects.append(f"v8_{key}_invalid")
    components = v8.get("momentum_component_relative")
    if (
        not isinstance(components, list)
        or len(components) != 3
        or not all(_positive_finite_number(x) for x in components)
    ):
        defects.append("v8_momentum_components_invalid")
    if not _positive_finite_number(last.get("mass_normalized_imbalance")):
        defects.append("mass_balance_missing_or_invalid")
    proof["checks"] = checks
    proof["defects"] = defects
    proof["nonlinear_v8_proven"] = not defects
    return proof
