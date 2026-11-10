"""Fail-closed, independently referenced *evidence* qualification for CFDX.

Registry PASS is historical metadata, not a verified result on this source.
One evidence JSON per registry ID is required, plus hashed raw solver artifacts.
This validates numerical *certificates*, not the correctness of the solver itself:
independent PDE/MMS and control-volume oracles are still required.

Use --require-all in a release qualification gate. With no evidence directory it
emits 25 explicit NOT_QUALIFIED records, never promotes registry status to PASS.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
from typing import Any

SCHEMA = "cfdx-physical-evidence/v1"
MAX_RESIDUAL = 1.0e-6  # upper bound, NOT a replacement for the case tolerance
# Mandatory independent observables for every registered qualification. A case
# with just one conveniently accurate scalar can no longer hide errors in the
# other claimed physics. References must supply values, units and uncertainties
# for ALL these obligations, not invented targets taken from solver output.
REQUIRED_QOI = {
    "LAM-COUETTE": ("u_mean", "velocity_profile_l2", "wall_shear"),
    "LAM-POISEUILLE": ("volume_flow", "velocity_profile_l2", "pressure_gradient"),
    "INC-GHIA": ("u_centerline_rms", "v_centerline_rms", "vortex_position_error"),
    "VER-MMS": ("field_l2", "field_linf", "face_flux_error"),
    "FORCE-VMFL036": ("Cd", "Cd_pressure", "Cd_viscous"),
    "EXT-NACA0012": ("Cd", "Cl", "Cm", "Cp_profile_l2", "Cf_profile_l2"),
    "INT-CHANNEL": ("bulk_velocity", "friction_factor", "pressure_gradient"),
    "SEP-BFS": ("reattachment_length", "velocity_profile_l2"),
    "TUR-FLATPLATE": ("Cf_profile_l2", "Uplus_profile_l2"),
    "TUR-CHANNEL": ("velocity_profile_l2", "reynolds_stress_l2"),
    "SEP-BFS-DRIVER": ("reattachment_length", "wall_shear_l2"),
    "TH-MIXED-CONVECTION": ("Nu_mean", "temperature_l2", "pressure_drop"),
    "TH-TJUNCTION": ("mixing_temperature_l2", "pressure_drop", "energy_balance"),
    "TUR-NACA0012": ("Cd", "Cl", "Cm", "Cp_profile_l2", "Cf_profile_l2"),
    "TH-CONV": ("Nu_mean", "temperature_l2"),
    "CHT-INTERFACE": ("interface_heat_flux", "interface_temperature_jump"),
    "RAD-P1": ("heat_flux", "incident_radiation"),
    "RAD-DOM": ("heat_flux", "angular_error"),
    "RAD-S2S": ("heat_flux", "view_factor_reciprocity_error"),
    "N9-S8-NACA0012": ("Cd", "Cl", "Cm", "Cp_profile_l2", "Cf_profile_l2"),
    "T1-SOLID-COND": ("heat_flux", "temperature_l2"),
    "T2-FLUID-ENERGY": ("energy_flux", "temperature_l2"),
    "INC-TRANSPOSE-STRESS-MMS": ("momentum_error_l2", "shear_traction_error_l2"),
    "SCALAR-TRANSPORT-NONORTH-AFFINE": ("scalar_error_l2", "flux_error_l2"),
    "WALL-DISTANCE-GEOMETRY": ("distance_l2", "distance_linf", "near_wall_linf"),
}
# Explicit physics obligations: no substring heuristics silently skip a balance
# on a case whose registry domain happens to use a different French/English label.
FLOW_CASES = {
    "LAM-COUETTE",
    "LAM-POISEUILLE",
    "INC-GHIA",
    "FORCE-VMFL036",
    "EXT-NACA0012",
    "INT-CHANNEL",
    "SEP-BFS",
    "TUR-FLATPLATE",
    "TUR-CHANNEL",
    "SEP-BFS-DRIVER",
    "TH-MIXED-CONVECTION",
    "TH-TJUNCTION",
    "TUR-NACA0012",
    "TH-CONV",
    "N9-S8-NACA0012",
    "T2-FLUID-ENERGY",
    "INC-TRANSPOSE-STRESS-MMS",
}
ENERGY_CASES = {
    "TH-MIXED-CONVECTION",
    "TH-TJUNCTION",
    "TH-CONV",
    "CHT-INTERFACE",
    "T1-SOLID-COND",
    "T2-FLUID-ENERGY",
    "RAD-P1",
    "RAD-DOM",
    "RAD-S2S",
}
RADIATION_CASES = {"RAD-P1", "RAD-DOM", "RAD-S2S"}
ALL_CHECKS = (
    "reference_artifact",
    "provenance",
    "raw_evidence",
    "case_contract",
    "solver",
    "residuals",
    "qoi_stability",
    "conservation",
    "force_decomposition",
    "reference",
    "mesh_series",
    "spatial_order",
)
FORCE_RE = ("FORCE-", "EXT-NACA", "TUR-NACA", "N9-S8-NACA")


def number(
    value: Any, name: str, *, positive: bool = False, nonnegative: bool = False
) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise TypeError(f"{name}: numeric value required")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name}: nonfinite value")
    if positive and result <= 0 or nonnegative and result < 0:
        raise ValueError(f"{name}: invalid sign")
    return result


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def verify_artifact(base: Path, artifact: Any, context: str) -> dict:
    if not isinstance(artifact, dict) or set(artifact) != {"path", "sha256"}:
        raise ValueError(f"{context}: path+sha256 required")
    rel = artifact["path"]
    checksum = artifact["sha256"]
    if not isinstance(rel, str) or not rel or Path(rel).is_absolute():
        raise ValueError(f"{context}: invalid relative artifact path")
    if (
        not isinstance(checksum, str)
        or len(checksum) != 64
        or any(c not in "0123456789abcdef" for c in checksum)
    ):
        raise ValueError(f"{context}: invalid sha256")
    root = base.resolve()
    candidate = (root / rel).resolve()
    if not candidate.is_relative_to(root) or not candidate.is_file():
        raise ValueError(f"{context}: missing/outside evidence directory")
    if candidate.stat().st_size == 0 or digest(candidate) != checksum:
        raise ValueError(f"{context}: empty or hash-mismatched artifact")
    return {"file": rel, "sha256": checksum, "bytes": candidate.stat().st_size}


def verify_provenance(evidence: dict, base: Path) -> dict:
    provenance = evidence["provenance"]
    if not isinstance(provenance, dict):
        raise TypeError("provenance: invalid structure")
    kind = provenance.get("source_kind")
    revision = provenance.get("source_revision")
    if kind == "git":
        if (
            not isinstance(revision, str)
            or len(revision) != 40
            or any(c not in "0123456789abcdef" for c in revision)
            or provenance.get("dirty") is not False
        ):
            raise ValueError("provenance: clean 40-digit git commit required")
    elif kind == "archive":
        if (
            not isinstance(revision, str)
            or len(revision) != 64
            or any(c not in "0123456789abcdef" for c in revision)
        ):
            raise ValueError("provenance: source archive SHA-256 required")
    else:
        raise ValueError("provenance: unknown source kind")
    artifacts = provenance.get("artifacts")
    if not isinstance(artifacts, dict):
        raise TypeError("provenance: missing artifacts")
    required = {"binary", "configuration", "solver_log", "reference", "oracle_source"}
    if kind == "archive":
        required.add("source_archive")
    if not required.issubset(artifacts):
        raise ValueError(
            f"provenance: missing artifacts {sorted(required - set(artifacts))}"
        )
    records = {
        name: verify_artifact(base, artifacts[name], name) for name in sorted(required)
    }
    if kind == "archive" and records["source_archive"]["sha256"] != revision:
        raise ValueError(
            "provenance: archive source_revision is not its actual SHA-256"
        )
    return records


def verify_raw_mesh_evidence(meshes: list[dict], base: Path, case_id: str) -> dict:
    """Reconcile the last accepted iterate with a separately stored CSV and oracle.

    A hash identifies the bytes but does NOT prove algorithm independence.
    Keep `oracle_source` SHA and numerical MMS in the release evidence.
    """
    results = []
    for mesh in meshes:
        name = mesh["id"]
        artifacts = mesh.get("artifacts")
        required = {"mesh", "history", "independent_oracle"}
        if not isinstance(artifacts, dict) or not required.issubset(artifacts):
            raise ValueError(
                f"raw_evidence: {name} missing mesh/history/oracle artifacts"
            )
        records = {
            key: verify_artifact(base, artifacts[key], f"{name}.{key}")
            for key in sorted(required)
        }
        histpath = (base / artifacts["history"]["path"]).resolve()
        with histpath.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            residuals = mesh["residuals"]
            required_csv = {"iteration", *residuals}
            if (
                not reader.fieldnames
                or len(set(reader.fieldnames)) != len(reader.fieldnames)
                or not required_csv.issubset(reader.fieldnames)
            ):
                raise ValueError(f"raw_evidence: {name} invalid history columns")
            previous, last = 0, None
            for row in reader:
                if None in row or any(v is None for v in row.values()):
                    raise ValueError(f"raw_evidence: {name} truncated CSV row")
                iteration = int(row["iteration"])
                if iteration <= previous:
                    raise ValueError(f"raw_evidence: {name} history not increasing")
                previous = iteration
                for key in residuals:
                    number(float(row[key]), f"{name}.history.{key}", nonnegative=True)
                last = row
            if last is None or previous != mesh["iterations"]:
                raise ValueError(
                    f"raw_evidence: {name} missing final accepted iteration"
                )
            for key, value in residuals.items():
                observed = float(last[key])
                expected = number(value, f"{name}.residuals.{key}", nonnegative=True)
                if not math.isclose(observed, expected, abs_tol=1e-13, rel_tol=1e-6):
                    raise ValueError(f"raw_evidence: {name}.{key} stale final residual")
        oraclepath = (base / artifacts["independent_oracle"]["path"]).resolve()
        oracle = json.loads(oraclepath.read_text(encoding="utf-8"))
        if (
            oracle.get("schema") != "cfdx-independent-fv-oracle/v1"
            or oracle.get("case_id") != case_id
            or oracle.get("mesh_id") != name
        ):
            raise ValueError(f"raw_evidence: {name} oracle case/mesh mismatch")
        for key, expected in mesh["independent_residuals"].items():
            observed = number(
                oracle.get("residuals", {}).get(key),
                f"{name}.oracle.{key}",
                nonnegative=True,
            )
            if not math.isclose(
                observed,
                number(expected, key, nonnegative=True),
                abs_tol=1e-13,
                rel_tol=1e-6,
            ):
                raise ValueError(
                    f"raw_evidence: {name} oracle residual mismatch: {key}"
                )
        if oracle.get("balances") != mesh["balances"]:
            raise ValueError(
                f"raw_evidence: {name} independent integral balances mismatch"
            )
        if "Cd_control_volume" in mesh.get("quantities", {}):
            value = number(
                oracle.get("Cd_control_volume"), f"{name}.oracle.Cd_control_volume"
            )
            if not math.isclose(
                value,
                mesh["quantities"]["Cd_control_volume"],
                abs_tol=1e-10,
                rel_tol=1e-8,
            ):
                raise ValueError(f"raw_evidence: {name} independent force mismatch")
        results.append({"mesh_id": name, "verified_artifacts": records})
    return {"raw_meshes": results}


def verify_contract(evidence: dict, reference: dict, case: dict) -> dict:
    if evidence.get("schema") != SCHEMA or evidence.get("case_id") != case["id"]:
        raise ValueError("case_contract: wrong schema or case ID")
    setup = evidence["setup"]
    if not isinstance(setup, dict) or not setup:
        raise ValueError("case_contract: missing physical setup")
    if setup.get("model") not in (
        "laminar",
        "RANS",
        "LES",
        "thermal",
        "CHT",
        "radiation",
        "multiphase",
        "scalar",
        "verification",
    ):
        raise ValueError("case_contract: unspecified physical model")
    for key in ("reynolds", "mach", "alpha_degrees"):
        if key in setup:
            number(setup[key], f"setup.{key}", nonnegative=(key != "alpha_degrees"))
    if "NACA" in case["id"]:
        for key in (
            "reynolds",
            "mach",
            "alpha_degrees",
            "density",
            "dynamic_viscosity",
            "reference_area",
            "reference_chord",
        ):
            number(
                setup.get(key),
                f"NACA setup.{key}",
                positive=key
                in (
                    "reynolds",
                    "density",
                    "dynamic_viscosity",
                    "reference_area",
                    "reference_chord",
                ),
            )
    if reference.get("case_id") != case["id"] or reference.get("setup") != setup:
        raise ValueError("case_contract: reference describes different physics")
    if (
        not isinstance(reference.get("citation"), str)
        or not reference["citation"].strip()
    ):
        raise ValueError("case_contract: independent reference citation required")
    if reference.get("independent_of_solver") is not True:
        raise ValueError("case_contract: no independent physical reference")
    return {"model": setup["model"], "reference": reference["citation"]}


def verify_solver(meshes: Any) -> dict:
    if not isinstance(meshes, list) or len(meshes) < 3:
        raise ValueError("solver: at least three meshes required")
    scheme = None
    seen = set()
    for mesh in meshes:
        if not isinstance(mesh, dict):
            raise TypeError("solver: malformed mesh record")
        name = mesh.get("id")
        if not isinstance(name, str) or not name or name in seen:
            raise ValueError("solver: missing/duplicate mesh ID")
        seen.add(name)
        if mesh.get("status") != "CONVERGED":
            raise ValueError(f"solver: {name} has no converged final state")
        if (
            not isinstance(mesh.get("iterations"), int)
            or isinstance(mesh["iterations"], bool)
            or mesh["iterations"] < 1
        ):
            raise ValueError(f"solver: {name} invalid accepted iteration count")
        number(mesh.get("runtime_seconds"), f"solver: {name}.runtime", positive=True)
        if not isinstance(mesh.get("final_scheme"), str) or not mesh["final_scheme"]:
            raise ValueError(f"solver: {name} missing final discretization")
        if scheme is None:
            scheme = mesh["final_scheme"]
        if scheme != mesh["final_scheme"]:
            raise ValueError("solver: discretization changes between meshes")
        if "NACA" in name.upper() or "NACA" in str(mesh.get("case", "")).upper():
            # Name-based checks are insufficient for this; case-ID enforcement below.
            pass
    return {"meshes": len(meshes), "final_scheme": scheme}


def verify_residuals(meshes: list[dict], case: dict) -> dict:
    case_id = case["id"]
    is_naca = "NACA" in case_id
    required: tuple[str, ...] = ()
    if case_id in FLOW_CASES:
        required += ("momentum", "continuity", "mass_global", "mass_local")
    if case_id in ENERGY_CASES:
        required += ("energy",)
    if case_id in RADIATION_CASES:
        required += ("radiation",)
    if case_id in {"SCALAR-TRANSPORT-NONORTH-AFFINE"}:
        required += ("scalar",)
    if case_id == "VER-MMS":
        required += ("field",)
    maxima = {key: 0.0 for key in required}
    for mesh in meshes:
        name = mesh["id"]
        residuals = mesh.get("residuals")
        limits = mesh.get("residual_limits")
        independent = mesh.get("independent_residuals")
        if not all(isinstance(v, dict) for v in (residuals, limits, independent)):
            raise ValueError(f"residuals: {name} missing independent residual sets")
        for key in required:
            value = number(
                residuals.get(key), f"{name}.residuals.{key}", nonnegative=True
            )
            independent_value = number(
                independent.get(key), f"{name}.independent.{key}", nonnegative=True
            )
            limit = number(limits.get(key), f"{name}.limit.{key}", positive=True)
            cap = 1e-8 if is_naca else MAX_RESIDUAL
            if limit > cap or value > limit or independent_value > limit:
                raise ValueError(
                    f"residuals: {name}.{key} exceeds strict verified limit"
                )
            maxima[key] = max(maxima[key], value, independent_value)
        if is_naca and (
            mesh["final_scheme"] not in (
                "TVD", "TVD_VANLEER", "TVD_VANALBADA", "MUSCL_CELL_LIMITED"
            )
            mesh["final_scheme"]
            not in ("TVD", "TVD_VANLEER", "TVD_VANALBADA", "MUSCL_CELL_LIMITED")
            or number(mesh.get("high_order_fraction"), "TVD fraction") != 1.0
        ):
            raise ValueError(f"residuals: {name} NACA must finish in full TVD")
    return {"max_independent_or_solver_residual": maxima}


def verify_qoi_stability(
    meshes: list[dict], reference: dict, case: dict, base: Path
) -> dict:
    """Recheck physical QoI drift from the actual hashed accepted-iterate CSV.

    Three recent accepted states are the minimum; use 10% of the independently
    declared reference error budget as the allowed instantaneous QoI drift.
    A low algebraic residual with a drifting Cd/Nu/profile never passes.
    """
    expected = reference.get("quantities", {})
    required = REQUIRED_QOI[case["id"]]
    checked = {}
    for mesh in meshes:
        name = mesh["id"]
        artifact = mesh.get("artifacts", {}).get("history", {})
        path = (base / artifact["path"]).resolve()
        if not path.is_relative_to(base.resolve()) or not path.is_file():
            raise ValueError(
                f"qoi_stability: {name} history outside evidence directory"
            )
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
                raise ValueError(f"qoi_stability: {name} lacks raw QoI history columns")
            last_three: list[dict] = []
            for row in reader:
                last_three.append(row)
                if len(last_three) > 3:
                    last_three.pop(0)
        if len(last_three) != 3:
            raise ValueError(f"qoi_stability: {name} requires >=3 accepted QoI samples")
        for key in required:
            spec = expected.get(key)
            if not isinstance(spec, dict):
                raise TypeError(
                    f"qoi_stability: {name} missing independent tolerance for {key}"
                )
            allowed = (
                number(spec.get("absolute_tolerance"), f"{key}.abs", nonnegative=True)
                + number(spec.get("relative_tolerance"), f"{key}.rel", nonnegative=True)
                * abs(number(spec.get("value"), f"{key}.reference"))
                + number(
                    spec.get("uncertainty"), f"{key}.uncertainty", nonnegative=True
                )
            )
            if allowed <= 0:
                raise ValueError(
                    f"qoi_stability: {key} needs positive reference error budget"
                )
            readings = [
                number(float(row[key]), f"{name}.history.{key}") for row in last_three
            ]
            final = number(mesh["quantities"].get(key), f"{name}.{key}")
            if not math.isclose(readings[-1], final, abs_tol=1e-12, rel_tol=1e-8):
                raise ValueError(f"qoi_stability: {name}.{key} final QoI is stale")
            drift = max(abs(value - readings[-1]) for value in readings[:-1])
            if drift > 0.1 * allowed:
                raise ValueError(
                    f"qoi_stability: {name}.{key} drift {drift:.4g} > {0.1 * allowed:.4g}"
                )
            checked[f"{name}.{key}"] = drift
    return {"max_absolute_drift_by_mesh_and_qoi": checked}


def verify_conservation(meshes: list[dict], case: dict) -> dict:
    case_id = case["id"]
    essential = set()
    if case_id in FLOW_CASES:
        essential.add("mass")
    if "NACA" in case_id:
        essential.update(("momentum_x", "momentum_y"))
    if case_id == "CHT-INTERFACE":
        essential.add("interface_energy")
    if case_id == "SCALAR-TRANSPORT-NONORTH-AFFINE":
        essential.add("scalar")
    if case_id in ENERGY_CASES:
        essential.add("energy")
    if case_id in RADIATION_CASES:
        essential.add("radiation")
    worst: dict[str, float] = {}
    for mesh in meshes:
        balances = mesh.get("balances")
        if not isinstance(balances, dict) or not essential.issubset(balances):
            raise ValueError(
                f"conservation: {mesh['id']} missing {sorted(essential)} balance"
            )
        for name, terms in balances.items():
            if not isinstance(terms, dict):
                raise TypeError(f"conservation: {name} bad balance")
            # Outward includes all boundary fluxes in the declared convention.
            values = {
                key: number(terms.get(key), f"{name}.{key}")
                for key in ("inward", "outward", "source", "storage")
            }
            limit = number(
                terms.get("relative_limit"), f"{name}.relative_limit", positive=True
            )
            if limit > MAX_RESIDUAL:
                raise ValueError(f"conservation: {name} tolerance too loose")
            scale = max(sum(abs(v) for v in values.values()), 1e-12)
            defect = (
                abs(
                    values["inward"]
                    + values["source"]
                    - values["outward"]
                    - values["storage"]
                )
                / scale
            )
            if defect > limit:
                raise ValueError(
                    f"conservation: {mesh['id']}.{name} defect {defect:.3g}"
                )
            worst[name] = max(worst.get(name, 0.0), defect)
    return {"max_relative_integral_defect": worst}


def verify_force(meshes: list[dict], case: dict) -> dict:
    if not case["id"].startswith(FORCE_RE):
        return {"applicable": False}
    for mesh in meshes:
        observable = mesh.get("quantities", {})
        if not isinstance(observable, dict):
            raise TypeError("force_decomposition: missing quantities")
        for key in ("Cd", "Cd_pressure", "Cd_viscous"):
            number(observable.get(key), f"{mesh['id']}.{key}")
        if not math.isclose(
            observable["Cd"],
            observable["Cd_pressure"] + observable["Cd_viscous"],
            abs_tol=1e-9,
            rel_tol=1e-7,
        ):
            raise ValueError(f"force_decomposition: {mesh['id']} Cd mismatch")
        # Independent integration evidence must be a distinct observable,
        # never the same solver accumulator reused twice.
        number(
            observable.get("Cd_independent_surface"),
            f"{mesh['id']}.Cd_independent_surface",
        )
        number(observable.get("Cd_control_volume"), f"{mesh['id']}.Cd_control_volume")
        for key in ("Cd_independent_surface", "Cd_control_volume"):
            if abs(observable["Cd"] - observable[key]) > 1e-3 * max(
                1.0, abs(observable["Cd"])
            ):
                raise ValueError(
                    f"force_decomposition: {mesh['id']} {key} disagreement"
                )
    return {"applicable": True, "force_checks": 3 * len(meshes)}


def verify_reference(meshes: list[dict], reference: dict, case: dict) -> dict:
    expected = reference.get("quantities")
    if not isinstance(expected, dict) or not expected:
        raise ValueError("reference: independent quantitative reference missing")
    required = set(REQUIRED_QOI[case["id"]])
    if not required.issubset(expected):
        raise ValueError(
            f"reference: missing independent QoIs {sorted(required - set(expected))}"
        )
    fine = min(meshes, key=lambda row: number(row.get("h"), "mesh.h", positive=True))
    errors = {}
    for key, spec in expected.items():
        if not isinstance(spec, dict):
            raise TypeError(f"reference: invalid {key} contract")
        if not isinstance(spec.get("unit"), str) or not spec["unit"].strip():
            raise ValueError(f"reference: missing {key} units")
        target = number(spec.get("value"), f"reference.{key}")
        uncertainty = number(
            spec.get("uncertainty"), f"reference.{key}.uncertainty", nonnegative=True
        )
        absolute = number(
            spec.get("absolute_tolerance"),
            f"reference.{key}.absolute_tolerance",
            nonnegative=True,
        )
        relative = number(
            spec.get("relative_tolerance"),
            f"reference.{key}.relative_tolerance",
            nonnegative=True,
        )
        if absolute == 0 and relative == 0 and uncertainty == 0:
            raise ValueError(f"reference: {key} zero tolerance/uncertainty")
        computed = number(fine.get("quantities", {}).get(key), f"{fine['id']}.{key}")
        error = abs(computed - target)
        allowed = absolute + relative * abs(target) + uncertainty
        if error > allowed:
            raise ValueError(
                f"reference: {key} physical error {error:.4g} > {allowed:.4g} {spec['unit']}"
            )
        errors[key] = {
            "computed": computed,
            "target": target,
            "error": error,
            "allowed": allowed,
            "unit": spec["unit"],
        }
    return {"fine_mesh": fine["id"], "quantities": errors}


def verify_mesh_series(meshes: list[dict], case: dict) -> dict:
    ordered = sorted(
        meshes,
        key=lambda item: number(item.get("h"), "mesh.h", positive=True),
        reverse=True,
    )
    hs = [number(row["h"], f"{row['id']}.h", positive=True) for row in ordered]
    if any(h0 <= h1 for h0, h1 in itertools.pairwise(hs)):
        raise ValueError("mesh_series: duplicate or nondecreasing spacing")
    refinement = [h0 / h1 for h0, h1 in itertools.pairwise(hs)]
    if any(r < 1.25 for r in refinement):
        raise ValueError("mesh_series: insufficient refinement ratio")
    if max(refinement) / min(refinement) > 1.05:
        raise ValueError("mesh_series: incompatible uniform Richardson ratios")
    family = {row.get("mesh_family") for row in ordered}
    if (
        None in family
        or len(family) != 1
        or not all(isinstance(name, str) and name.strip() for name in family)
    ):
        raise ValueError("mesh_series: unverified or mixed mesh families")
    wall_layers = []
    max_nonorth = []
    for row in ordered:
        quality = row.get("mesh_quality")
        if not isinstance(quality, dict):
            raise TypeError("mesh_series: missing mesh quality")
        if (
            quality.get("inverted_cells") != 0
            or number(
                quality.get("minimum_cell_volume"),
                "mesh_quality.minimum_cell_volume",
                positive=True,
            )
            <= 0
        ):
            raise ValueError("mesh_series: inverted or degenerate cells")
        max_angle = number(
            quality.get("max_nonorthogonality_degrees"),
            "mesh_quality.max_nonorthogonality",
            nonnegative=True,
        )
        if max_angle >= 80:
            raise ValueError("mesh_series: near-degenerate nonorthogonality")
        max_nonorth.append(max_angle)
        if "NACA" in case["id"]:
            for key in (
                "minimum_first_layer_spacing",
                "mean_first_layer_spacing",
                "maximum_first_layer_spacing",
                "p95_nonorthogonality_degrees",
                "maximum_cell_edge_aspect_ratio",
            ):
                number(
                    quality.get(key),
                    f"{row['id']}.{key}",
                    positive=key != "p95_nonorthogonality_degrees",
                    nonnegative=key == "p95_nonorthogonality_degrees",
                )
            wall_layers.append(quality)
    if "NACA" in case["id"]:
        for i, (r, a, b) in enumerate(zip(refinement, wall_layers, wall_layers[1:])):
            for key in (
                "minimum_first_layer_spacing",
                "mean_first_layer_spacing",
                "maximum_first_layer_spacing",
            ):
                actual_ratio = a[key] / b[key]
                if abs(actual_ratio / r - 1) > 0.2:
                    raise ValueError(
                        f"mesh_series: {key} does not follow the tangential refinement"
                    )
            if (
                abs(max_nonorth[i + 1] - max_nonorth[i]) > 10
                or abs(
                    a["p95_nonorthogonality_degrees"]
                    - b["p95_nonorthogonality_degrees"]
                )
                > 10
            ):
                raise ValueError(
                    "mesh_series: NACA nonorthogonality distribution drifts"
                )
    return {
        "order": [row["id"] for row in ordered],
        "h": hs,
        "ratios": refinement,
        "family": next(iter(family)),
        "max_nonorthogonality_degrees": max_nonorth,
    }


def verify_order(meshes: list[dict], reference: dict) -> dict:
    ordered = sorted(meshes, key=lambda row: float(row["h"]), reverse=True)
    if len(ordered) != 3:
        raise ValueError("spatial_order: use exactly three meshes for verified order")
    h = [float(row["h"]) for row in ordered]
    ratio = h[1] / h[2]
    detail = {}
    for key in reference["quantities"]:
        f = [number(row["quantities"].get(key), f"{key}.grid") for row in ordered]
        d01, d12 = f[0] - f[1], f[1] - f[2]
        if d01 == 0 and d12 == 0:
            # Identical values on three meshes cannot establish observed order;
            # an exact analytic solution needs a SEPARATE MMS/round-off proof.
            raise ValueError(
                f"spatial_order: {key} identical-grid plateau has no observed GCI"
            )
        if d01 * d12 <= 0 or abs(d12) >= abs(d01):
            raise ValueError(f"spatial_order: {key} oscillating/stagnating mesh series")
        r21, r32 = h[0] / h[1], h[1] / h[2]
        log21, log32 = math.log(r21), math.log(r32)
        factor = abs(d01 / d12)

        def diff_ratio(p, log21=log21, log32=log32):
            return math.expm1(p * log21) / (-math.expm1(-p * log32))

        if factor <= log21 / log32 or factor >= diff_ratio(8):
            raise ValueError(f"spatial_order: {key} invalid positive observed order")
        lo, hi = 1e-9, 8.0
        for _ in range(70):
            mid = (lo + hi) / 2
            if diff_ratio(mid) < factor:
                lo = mid
            else:
                hi = mid
        p = (lo + hi) / 2
        # Require an appreciable order; a 0.001 order 'GCI' is not credible.
        if p < 0.5:
            raise ValueError(f"spatial_order: {key} order {p:.3g} < 0.5")
        extrapolated = f[2] + d12 * (-1) / math.expm1(p * math.log(ratio))
        spec = reference["quantities"][key]
        target = number(spec["value"], f"reference.{key}")
        allowed = (
            number(spec["absolute_tolerance"], "absolute_tolerance", nonnegative=True)
            + number(spec["relative_tolerance"], "relative_tolerance", nonnegative=True)
            * abs(target)
            + number(spec["uncertainty"], "uncertainty", nonnegative=True)
        )
        gci_absolute = 1.25 * abs(d12) / math.expm1(p * math.log(ratio))
        if gci_absolute > allowed:
            raise ValueError(
                f"spatial_order: {key} fine-grid GCI {gci_absolute:.4g} > physical tolerance {allowed:.4g}"
            )
        if abs(extrapolated - target) > allowed:
            raise ValueError(
                f"spatial_order: {key} extrapolated value is outside physical reference"
            )
        detail[key] = {
            "status": "RICHARDSON_CANDIDATE",
            "p": p,
            "extrapolated": extrapolated,
            "fine_grid_absolute_gci": gci_absolute,
        }
    return detail


def evaluate_case(case: dict, path: Path) -> dict:
    result = {
        "id": case["id"],
        "registry_status": case.get("status"),
        "verification_status": "NOT_QUALIFIED",
        "checks": [],
    }
    if not path.is_file():
        result["checks"].append(
            {
                "check": "evidence_file",
                "status": "MISSING",
                "reason": "No solver/references file on this HEAD",
            }
        )
        result["checks"].extend(
            {"check": name, "status": "MISSING"} for name in ALL_CHECKS
        )
        result["required_physical_observables"] = list(REQUIRED_QOI[case["id"]])
        return result
    try:
        evidence = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError, TypeError) as error:
        result["checks"].append(
            {
                "check": "evidence_file",
                "status": "FAIL",
                "reason": f"Unreadable evidence: {error}",
            }
        )
        return result
    if not isinstance(evidence, dict):
        result["checks"].append(
            {
                "check": "evidence_file",
                "status": "FAIL",
                "reason": "Evidence must be an object",
            }
        )
        return result

    def gate(name, action):
        try:
            details = action()
            result["checks"].append(
                {"check": name, "status": "PASS", "details": details}
            )
        except (ValueError, KeyError, OSError, TypeError, OverflowError) as error:
            result["checks"].append(
                {"check": name, "status": "FAIL", "reason": str(error)}
            )

    base = path.parent
    meshes = evidence.get("meshes", [])
    # Independent artifact: exact bytes hash-matched, not a self-referential
    # reference JSON encoded in the same solver-generated summary.
    reference_path = base / "reference.json"

    def read_reference():
        if not reference_path.is_file():
            raise ValueError("missing independent reference.json")
        ref = json.loads(reference_path.read_text(encoding="utf-8"))
        expected_sha = evidence["provenance"]["artifacts"]["reference"]["sha256"]
        if digest(reference_path) != expected_sha:
            raise ValueError("reference SHA mismatch")
        return ref

    try:
        reference = read_reference()
    except (ValueError, KeyError, OSError, TypeError) as error:
        reference = None
        result["checks"].append(
            {"check": "reference_artifact", "status": "FAIL", "reason": str(error)}
        )
    else:
        result["checks"].append({"check": "reference_artifact", "status": "PASS"})
    gate("provenance", lambda: verify_provenance(evidence, base))
    gate("raw_evidence", lambda: verify_raw_mesh_evidence(meshes, base, case["id"]))
    gate("case_contract", lambda: verify_contract(evidence, reference, case))
    gate("solver", lambda: verify_solver(meshes))
    gate("residuals", lambda: verify_residuals(meshes, case))
    gate("qoi_stability", lambda: verify_qoi_stability(meshes, reference, case, base))
    gate("conservation", lambda: verify_conservation(meshes, case))
    gate("force_decomposition", lambda: verify_force(meshes, case))
    gate("reference", lambda: verify_reference(meshes, reference, case))
    gate("mesh_series", lambda: verify_mesh_series(meshes, case))
    gate("spatial_order", lambda: verify_order(meshes, reference))
    if all(row["status"] == "PASS" for row in result["checks"]):
        result["verification_status"] = "EVIDENCE_GATES_PASS"
    return result


def evaluate_registry(
    registry: dict, evidence_dir: Path, selection: set[str] | None = None
) -> dict:
    cases = registry.get("cases", [])
    if not isinstance(cases, list) or not cases:
        raise ValueError("empty qualification registry")
    ids = [case.get("id") for case in cases]
    if not all(isinstance(i, str) and i for i in ids) or len(set(ids)) != len(ids):
        raise ValueError("invalid or duplicate qualification IDs")
    if set(ids) != set(REQUIRED_QOI):
        raise ValueError(
            f"equation/physical obligation map does not cover registry IDs: "
            f"missing={sorted(set(ids) - set(REQUIRED_QOI))}, "
            f"retired={sorted(set(REQUIRED_QOI) - set(ids))}"
        )
    if selection is not None and not selection.issubset(ids):
        raise ValueError(f"unknown case IDs: {sorted(selection - set(ids))}")
    selected = [case for case in cases if selection is None or case["id"] in selection]
    # Do not use unchecked IDs as paths; case IDs come from a JSON file.
    for case in selected:
        if not all(c.isascii() and (c.isalnum() or c in "_-") for c in case["id"]):
            raise ValueError(f"unsafe case ID: {case['id']!r}")
    reports = [
        evaluate_case(case, evidence_dir / case["id"] / "evidence.json")
        for case in selected
    ]
    passed = sum(r["verification_status"] == "EVIDENCE_GATES_PASS" for r in reports)
    return {
        "schema": "cfdx-qualification-coverage/v1",
        "cases": reports,
        "coverage": {
            "selected": len(reports),
            "evidence_pass": passed,
            "not_qualified": len(reports) - passed,
        },
        "all_evidence_gates_pass": passed == len(reports),
        "warning": "Evidence gates are not proof of equation correctness: independent MMS and experiments remain necessary",
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--registry",
        type=Path,
        default=root / "docs/validation/CFDX_QUALIFICATION_REGISTRY.json",
    )
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--case", action="append", dest="cases")
    parser.add_argument(
        "--require-all", action="store_true", help="exit 2 if any case unqualified"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        registry = json.loads(args.registry.read_text(encoding="utf-8"))
        report = evaluate_registry(
            registry, args.evidence_dir, set(args.cases) if args.cases else None
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temp = args.output.with_name(args.output.name + ".tmp")
        temp.write_text(
            json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temp, args.output)
    except (OSError, ValueError, KeyError, TypeError) as error:
        args.output.unlink(missing_ok=True)
        print(f"QUALIFICATION_EVIDENCE_GATE: FAIL {error}")
        return 2
    print(
        f"QUALIFICATION_EVIDENCE_GATE: {report['coverage']['evidence_pass']}/"
        f"{report['coverage']['selected']} verifiable; report={args.output}"
    )
    return 2 if args.require_all and not report["all_evidence_gates_pass"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
