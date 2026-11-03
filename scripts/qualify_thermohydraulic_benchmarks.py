"""Fail-closed postprocessor for independently referenced TH-01/03/07/08 CFDX runs.

This program deliberately DOES NOT run the solver. The input must contain actual
production CFDX outputs at three mesh levels, with independently sourced,
SHA-256-pinned reference profiles. With no solver evidence the physical status
is INCOMPLETE, never PASS. See the benchmark README for the input contract.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
CONTRACTS = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
CASES = ("TH-01", "TH-03", "TH-07", "TH-08")
SCENARIOS = {
    "TH-01": ("Re_tau_180", "Re_tau_395", "Re_tau_590"),
    "TH-03": ("Re_h_37500",),
    "TH-07": ("Gr_0", "Gr_640000", "Gr_950000", "Gr_1600000"),
    "TH-08": ("isothermal", "non_isothermal"),
}
# Quantities required in a steady mean-field qualification. These are NOT the
# complete DNS/LES fluctuation gates; only a time-resolved run can claim those.
PROFILE_NAMES = {
    "TH-01": {s: ("u_plus", "reynolds_uv_plus", "total_shear_plus") for s in SCENARIOS["TH-01"]},
    "TH-03": {"Re_h_37500": ("cf", "cp", "u_mean")},
    "TH-07": {s: ("u_mean", "temperature_mean") for s in SCENARIOS["TH-07"]},
    "TH-08": {
        "isothermal": ("u_mean",),
        "non_isothermal": ("u_mean", "temperature_mean"),
    },
}
# Predeclared numerical limits, NOT post-fit to CFDX data. These only establish
# readiness for meaningful physical comparison, not model-form validation.
MAX_MASS_IMBALANCE = 1.0e-4
MAX_ENERGY_IMBALANCE = 1.0e-3
MAX_TRUE_RESIDUAL = 1.0e-6
MAX_PROFILE_WEIGHTED_RMS = 2.0
MAX_PROFILE_WEIGHTED_PEAK = 5.0
MAX_MESH_QOI_DRIFT = 0.05
MAX_PROFILE_MESH_DRIFT = 0.05
# Canonical exported normalizations. A mixing of dimensional and wall units
# can produce an apparently excellent numerical fit with meaningless values.
PROFILE_UNITS = {
    "TH-01": {"u_plus": "U/u_tau", "reynolds_uv_plus": "-uv/u_tau^2",
              "total_shear_plus": "tau/tau_wall"},
    "TH-03": {"cf": "Cf", "cp": "Cp", "u_mean": "U/Uref"},
    "TH-07": {"u_mean": "U/Ub", "temperature_mean": "Theta"},
    "TH-08": {"u_mean": "U/Ub", "temperature_mean": "Tstar"},
}


class EvidenceError(ValueError):
    """Missing, malformed, or physically inadmissible qualification evidence."""


def check(ok: bool, message: str) -> None:
    if not ok:
        raise EvidenceError(message)


def number(value: Any, name: str, *, nonnegative: bool = False) -> float:
    check(type(value) in (int, float), f"{name}: expected numeric value")
    result = float(value)
    check(math.isfinite(result), f"{name}: nonfinite value")
    if nonnegative:
        check(result >= 0.0, f"{name}: negative value")
    return result


def checked_file(base: Path, relative: str, sha256: str) -> Path:
    check(isinstance(relative, str) and relative != "", "file path missing")
    check(isinstance(sha256, str) and re.fullmatch(r"[0-9a-f]{64}", sha256) is not None,
          f"{relative}: SHA-256 must be pinned in lowercase hex")
    path = (base / relative).resolve()
    check(path.is_relative_to(base.resolve()), f"{relative}: file outside evidence directory")
    check(path.is_file(), f"{relative}: missing file")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    check(digest == sha256, f"{relative}: SHA-256 mismatch (reference changed)")
    return path


def read_points(path: Path, *, reference: bool) -> list[tuple[float, float, float]]:
    """Read x,value[,sigma] as comma-separated data without extrapolation."""
    rows: list[tuple[float, float, float]] = []
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        needed = {"x", "value", "sigma"} if reference else {"x", "value"}
        check(reader.fieldnames is not None and needed.issubset(reader.fieldnames),
              f"{path.name}: expected CSV columns {sorted(needed)}")
        for ix, row in enumerate(reader, 2):
            try:
                x = float(row["x"])
                y = float(row["value"])
                sigma = float(row["sigma"]) if reference else 0.0
            except (ValueError, TypeError, KeyError) as exc:
                raise EvidenceError(f"{path.name}:{ix}: malformed CSV row") from exc
            check(all(map(math.isfinite, (x, y, sigma))), f"{path.name}:{ix}: nonfinite")
            if reference:
                check(sigma > 0.0, f"{path.name}:{ix}: sigma must be positive")
            if rows:
                check(x > rows[-1][0], f"{path.name}:{ix}: x must increase strictly")
            rows.append((x, y, sigma))
    check(len(rows) >= 5, f"{path.name}: at least five independent profile points required")
    return rows


def interpolate(points: list[tuple[float, float, float]], x: float) -> float:
    check(points[0][0] <= x <= points[-1][0], "profile extrapolation forbidden")
    for i, (xa, ya, _) in enumerate(points[:-1]):
        xb, yb, _ = points[i + 1]
        if xa <= x <= xb:
            return ya + (yb - ya) * (x - xa) / (xb - xa)
    return points[-1][1]


def compare_profile(base: Path, case: str, entry: dict[str, Any]) -> dict[str, Any]:
    name = entry.get("name")
    check(isinstance(name, str), "profile name missing")
    check(bool(entry.get("reference_source_url")), f"{name}: reference provenance missing")
    check(bool(entry.get("station")) and bool(entry.get("units")),
          f"{name}: explicit station and physical units/normalization required")
    check(entry["units"] == PROFILE_UNITS[case].get(name),
          f"{name}: unsupported normalization {entry['units']!r} for {case}")
    reference = read_points(checked_file(base, entry.get("reference_file", ""),
                                         entry.get("reference_sha256", "")), reference=True)
    simulation = read_points(checked_file(base, entry.get("simulation_file", ""),
                                          entry.get("simulation_sha256", "")), reference=False)
    # Require every published station to lie in the CFDX coordinate interval.
    check(simulation[0][0] <= reference[0][0] and
          simulation[-1][0] >= reference[-1][0],
          f"{name}: CFDX does not cover full reference profile (no extrapolation)")
    sampled = [interpolate(simulation, x) for x, _, _ in reference]
    errors = [(predicted - y) / sigma for predicted, (_, y, sigma) in zip(sampled, reference, strict=True)]
    wrms = math.sqrt(sum(e * e for e in errors) / len(errors))
    wpeak = max(abs(e) for e in errors)
    check(wrms <= MAX_PROFILE_WEIGHTED_RMS and wpeak <= MAX_PROFILE_WEIGHTED_PEAK,
          f"{name}: physical mismatch: weighted RMS={wrms:.4g}, peak={wpeak:.4g}")
    return {"name": name, "station": entry["station"], "units": entry["units"],
            "n_reference_points": len(errors), "weighted_rms": wrms,
            "weighted_peak": wpeak, "comparison_status": "WITHIN_FIXED_LIMITS",
            "_sampled_simulation": sampled, "_reference_sha256": entry["reference_sha256"]}


def scalar_targets(case: str, scenario: str) -> dict[str, tuple[float, float]]:
    data = json.loads((CONTRACTS / f"{case}.json").read_text(encoding="utf-8"))
    if case == "TH-01":
        nominal = int(scenario.removeprefix("Re_tau_"))
        row = next(x for x in data["physics"]["reference_cases"] if x["re_tau_nominal"] == nominal)
        # This is a flow-state matching bound, NOT an accuracy score for U+.
        return {"re_tau": (row["re_tau_actual"], 0.03 * row["re_tau_actual"])}
    if case == "TH-03":
        qoi = data["qoi"]["reference_targets"]
        return {"reattachment_x_over_h": (qoi["reattachment_x_over_h"],
                                            qoi["reattachment_uncertainty_absolute"])}
    if case == "TH-07":
        row = {"Gr_0": "Gr_0", "Gr_640000": "Gr_6.4e5",
               "Gr_950000": "Gr_9.5e5", "Gr_1600000": "Gr_1.6e6"}[scenario]
        stress = data["reference_targets"]["wall_shear_stress_mean"][row]
        if scenario == "Gr_0":
            return {"tau_wall": (stress["reference_isothermal"], 0.05 * stress["reference_isothermal"])}
        return {wall: (stress[wall], 0.05 * stress[wall]) for wall in ("heated", "cooled")}
    if case == "TH-08" and scenario == "non_isothermal":
        return {"t_star_mixed": (data["physics"]["derived"]["mixed_temperature_star_fully_mixed"],
                                  data["measurement_uncertainty"]["normalized_mean_temperature_absolute"])}
    return {}


def check_native_run_evidence(base: Path, case: str, scenario: str,
                              run: dict[str, Any], tag: str) -> None:
    """Verify a native production trace, not the declarations in a profile JSON.

    The experimental runner writes run_evidence.json next to its unmodified
    CFDX outputs and records their SHA-256. A manually assembled profile
    manifest is not sufficient to satisfy the strict production gate.
    """
    evidence_path = checked_file(base, run.get("native_run_evidence_file", ""),
                                 run.get("native_run_evidence_sha256", ""))
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    check(evidence.get("status") == "EXECUTED_UNQUALIFIED" and
          evidence.get("physical_qualification") == "NOT_QUALIFIED" and
          evidence.get("solver_executed") is True and
          type(evidence.get("exit_code")) is int and evidence["exit_code"] == 0,
          f"{tag}: native production execution failed or missing")
    check((evidence.get("case_id"), evidence.get("scenario_id"),
           evidence.get("mesh_id")) == (case, scenario, run["mesh_id"]),
          f"{tag}: native run case/scenario/mesh mismatch")
    check(isinstance(evidence.get("solver_binary_sha256"), str) and
          re.fullmatch(r"[0-9a-f]{64}", evidence["solver_binary_sha256"]) is not None and
          isinstance(evidence.get("case_file_sha256"), str) and
          re.fullmatch(r"[0-9a-f]{64}", evidence["case_file_sha256"]) is not None,
          f"{tag}: native binary or mesh hash absent")
    checked_file(evidence_path.parent, "solver.log", evidence.get("solver_log_sha256", ""))
    file_hashes = evidence.get("files_sha256")
    check(isinstance(file_hashes, dict), f"{tag}: native artifact checksums absent")
    expected_files = {
        "execution.json", "convergence.json",
        "benchmark_observables/mass_conservation.json",
        "benchmark_observables/cell_fields.csv",
        "benchmark_observables/face_fluxes.csv",
    }
    if case == 'TH-07':
        expected_files.update(('energy_convergence.json', 'energy_summary.json'))
    if case == 'TH-03':
        expected_files.add('wall_resolution.json')
    check(expected_files.issubset(file_hashes),
          f"{tag}: missing native solver-owned outputs")
    for name in expected_files:
        checked_file(evidence_path.parent, name, file_hashes[name])
    native_execution = json.loads((evidence_path.parent / "execution.json").read_text(encoding="utf-8"))
    native_convergence = json.loads((evidence_path.parent / "convergence.json").read_text(encoding="utf-8"))
    native_mass = json.loads((evidence_path.parent / "benchmark_observables/mass_conservation.json").read_text(encoding="utf-8"))
    check(native_execution.get("format") == "CFDX-EXECUTION" and
          native_execution.get("process_exit_code") == 0 and
          native_execution.get("converged") is True and
          native_convergence.get("format") == "CFDX-CONVERGENCE" and
          native_convergence.get("converged") is True and
          isinstance(native_convergence.get("history"), list) and native_convergence["history"] and
          native_mass.get("format") == "CFDX-THERMOHYDRAULIC-CONSERVATION" and
          native_mass.get("physical_qualification") == "NOT_QUALIFIED" and
          native_mass.get("cells") == run["cells"],
          f"{tag}: native solver status or cell count mismatch")
    last = native_convergence["history"][-1]
    true_residual = max(
        number(last.get("momentum_equation_residual_rhs_relative"),
               f"{tag}.native_momentum_residual", nonnegative=True),
        number(last.get("continuity_normalized"),
               f"{tag}.native_continuity_residual", nonnegative=True))
    mass = number(native_mass.get("global_mass_imbalance_relative"),
                  f"{tag}.native_mass_imbalance", nonnegative=True)
    max_cell_mass = number(native_mass.get("maximum_cell_mass_imbalance_relative"),
                           f"{tag}.native_max_cell_mass_imbalance", nonnegative=True)
    check(max_cell_mass <= MAX_MASS_IMBALANCE,
          f"{tag}: local cell mass leakage {max_cell_mass} exceeds {MAX_MASS_IMBALANCE}")
    if case == 'TH-03':
        try:
            from scripts.audit_th03_wall_resolution import audit as audit_sst_wall
        except ImportError:
            from audit_th03_wall_resolution import audit as audit_sst_wall
        native_wall = json.loads((evidence_path.parent/'wall_resolution.json').read_text())
        rho = number(native_wall.get('density_kg_m3'), f'{tag}.density')
        nu = number(native_wall.get('molecular_kinematic_viscosity_m2_s'),
                    f'{tag}.molecular_nu')
        check(rho > 0 and nu > 0, f'{tag}: invalid wall-audit physical scales')
        walls = native_wall.get('audited_wall_patches')
        check(isinstance(walls,list) and walls and all(isinstance(w,str) and w for w in walls),
              f'{tag}: missing all-wall coverage')
        derived = audit_sst_wall(
            evidence_path.parent/'benchmark_observables/face_fluxes.csv',
            evidence_path.parent/'benchmark_observables/cell_fields.csv',
            wall_patches=tuple(walls),density=rho,viscosity=nu)
        check(native_wall.get('status') == 'PREFLIGHT_PASS' and
              derived['status'] == 'PREFLIGHT_PASS' and
              native_wall.get('physical_qualification') == 'NOT_QUALIFIED' and
              native_wall.get('wall_method') == derived['wall_method'] and
              native_wall.get('n_faces') == derived['n_faces'] and
              abs(number(native_wall.get('max_y_plus_estimate'),f'{tag}.wall_y_plus') -
                  derived['max_y_plus_estimate']) < 1e-9,
              f'{tag}: native SST wall-resolution diagnostic does not pass')
    if case == 'TH-07':
        energy = json.loads((evidence_path.parent / 'energy_convergence.json').read_text())
        summary = json.loads((evidence_path.parent / 'energy_summary.json').read_text())
        check(energy.get('converged') is True and
              energy.get('linear_solver_converged') is True and
              energy.get('global_energy_balance_passed') is True and
              energy.get('local_energy_balance_passed') is True and
              energy.get('temperature_admissibility_passed') is True,
              f'{tag}: native thermal convergence contract incomplete')
        global_energy = number(energy.get('global_balance_normalized'),
                               f'{tag}.native_energy',nonnegative=True)
        local_energy = number(energy.get('max_cell_residual_normalized'),
                              f'{tag}.native_energy_cell',nonnegative=True)
        check(max(global_energy,local_energy) <= MAX_ENERGY_IMBALANCE,
              f'{tag}: native energy conservation exceeds gate')
        tmin = number(summary.get('min_temperature'),f'{tag}.temperature_min_K')
        tmax = number(summary.get('max_temperature'),f'{tag}.temperature_max_K')
        check(0<tmin<=tmax, f'{tag}: invalid native temperature envelope')
        declared = number(run.get('energy_imbalance_relative'),
                          f'{tag}.declared_energy',nonnegative=True)
        check(abs(declared-global_energy)<=1e-12*max(1,abs(global_energy)) and
              number(evidence.get('energy_imbalance_relative'),
                     f'{tag}.native_manifest_energy',nonnegative=True)==global_energy,
              f'{tag}: declared energy balance does not match solver-owned audit')
    for name, actual in (("true_residual_max", true_residual),
                         ("mass_imbalance_relative", mass)):
        declared = number(run.get(name), f"{tag}.{name}", nonnegative=True)
        check(abs(declared - actual) <= 1.e-12 * max(1., abs(actual)),
              f"{tag}: {name} is not backed by the native CFDX outputs")
    check(number(evidence.get("true_residual_max"), f"{tag}.native_run_residual") ==
          true_residual and
          number(evidence.get("mass_imbalance_relative"), f"{tag}.native_run_mass") == mass,
          f"{tag}: native run manifest contradicts saved solver output")


def check_meshes(case: str, scenario: str, runs: list[Any], base: Path,
                 *, require_native_evidence: bool = False) -> dict[str, Any]:
    check(isinstance(runs, list) and len(runs) >= 3,
          f"{scenario}: at least three distinct mesh levels required")
    targets = scalar_targets(case, scenario)
    meshes: list[dict[str, Any]] = []
    measurements: dict[str, list[float]] = {q: [] for q in targets}
    mesh_profiles: list[dict[tuple[str, str], tuple[str, list[float]]]] = []
    for i, run in enumerate(runs):
        check(isinstance(run, dict), f"{scenario}: malformed run")
        tag = f"{scenario}/mesh[{i}]"
        check(run.get("solver_executed") is True and run.get("converged") is True,
              f"{tag}: missing actual converged CFDX production run")
        check(isinstance(run.get("mesh_id"), str) and bool(run["mesh_id"]),
              f"{tag}: mesh ID missing")
        if require_native_evidence:
            check_native_run_evidence(base, case, scenario, run, tag)
        cells = number(run.get("cells"), f"{tag}.cells")
        h = number(run.get("characteristic_h"), f"{tag}.characteristic_h")
        check(cells >= 100 and h > 0, f"{tag}: invalid mesh resolution")
        if meshes:
            check(cells > meshes[-1]["cells"] and h < meshes[-1]["characteristic_h"],
                  f"{tag}: mesh ordering/refinement invalid")
            check(h <= 0.90 * meshes[-1]["characteristic_h"],
                  f"{tag}: insufficient systematic refinement")
        residual = number(run.get("true_residual_max"), f"{tag}.true_residual_max", nonnegative=True)
        mass = number(run.get("mass_imbalance_relative"), f"{tag}.mass_imbalance_relative", nonnegative=True)
        check(residual <= MAX_TRUE_RESIDUAL and mass <= MAX_MASS_IMBALANCE,
              f"{tag}: unconverged or mass imbalance too high: residual={residual}, mass={mass}")
        if case in ("TH-07", "TH-08") and scenario != "isothermal":
            energy = number(run.get("energy_imbalance_relative"),
                            f"{tag}.energy_imbalance_relative", nonnegative=True)
            check(energy <= MAX_ENERGY_IMBALANCE,
                  f"{tag}: energy imbalance {energy} > {MAX_ENERGY_IMBALANCE}")
        qoi = run.get("qoi")
        check(isinstance(qoi, dict), f"{tag}: missing measured QoIs")
        for name, (target, tolerance) in targets.items():
            actual = number(qoi.get(name), f"{tag}.{name}")
            measurements[name].append(actual)
            if i == len(runs) - 1:
                check(abs(actual - target) <= tolerance,
                      f"{tag}: {name}={actual} outside reference {target} +/- {tolerance}")
        profiles = run.get("profiles")
        check(isinstance(profiles, list), f"{tag}: missing profile comparisons")
        required = set(PROFILE_NAMES[case][scenario])
        present = {p.get("name") for p in profiles if isinstance(p, dict)}
        check(required.issubset(present),
              f"{tag}: missing mandatory profiles {sorted(required - present)}")
        station_keys = {(p.get("name"), p.get("station")) for p in profiles
                        if isinstance(p, dict)}
        check(len(station_keys) == len(profiles),
              f"{tag}: duplicate/malformed profile stations")
        if case == "TH-03":
            # Five published NASA/Driver-Seegmiller velocity stations;
            # profiles are station-resolved, not a single averaged curve.
            required_u_stations = {"x/h=-4", "x/h=1", "x/h=4", "x/h=6", "x/h=10"}
            seen_u = {p.get("station") for p in profiles if p.get("name") == "u_mean"}
            check(required_u_stations.issubset(seen_u),
                  f"{tag}: missing BFS velocity stations "
                  f"{sorted(required_u_stations - seen_u)}")
        comparisons = [compare_profile(base, case, profile) for profile in profiles]
        normalized_profiles = {}
        for comparison in comparisons:
            key = (comparison["name"], comparison["station"])
            normalized_profiles[key] = (comparison.pop("_reference_sha256"),
                                        comparison.pop("_sampled_simulation"))
        if mesh_profiles:
            check(set(normalized_profiles) == set(mesh_profiles[0]),
                  f"{tag}: profile stations change with mesh refinement")
            for key, (digest, samples) in normalized_profiles.items():
                original_digest, original_samples = mesh_profiles[0][key]
                check(digest == original_digest and len(samples) == len(original_samples),
                      f"{tag}: reference dataset or sampling positions change across meshes: {key}")
        mesh_profiles.append(normalized_profiles)
        meshes.append({"mesh_id": run["mesh_id"], "cells": int(cells),
                       "characteristic_h": h, "true_residual_max": residual,
                       "mass_imbalance_relative": mass, "profiles": comparisons})
    # A single grid-converged integral cannot certify a wrong or oscillating
    # velocity/temperature profile. Enforce mesh stability for EVERY required
    # station, including scenarios such as TH-08 isothermal with no scalar QoI.
    profile_stability = {}
    for key in sorted(mesh_profiles[-1]):
        a, b, c = (level[key][1] for level in mesh_profiles[-3:])
        scale = max(max(abs(v) for v in c), 1.0e-12)
        delta_coarse = math.sqrt(sum((v-w)**2 for v, w in zip(a, b, strict=True)) / len(c))
        delta_fine = math.sqrt(sum((v-w)**2 for v, w in zip(b, c, strict=True)) / len(c))
        drift = delta_fine / scale
        check(drift <= MAX_PROFILE_MESH_DRIFT,
              f"{scenario}: profile {key} still mesh-sensitive ({drift:.3%} > 5%)")
        check(delta_fine <= 1.0e-12 * scale or delta_coarse > delta_fine,
              f"{scenario}: profile {key} does not improve with refinement")
        profile_stability[f"{key[0]}@{key[1]}"] = {
            "fine_medium_relative_l2": drift,
            "coarse_medium_relative_l2": delta_coarse / scale,
        }
    stability: dict[str, float] = {}
    observed: dict[str, Any] = {}
    for name, values in measurements.items():
        # Refinement requires the last two levels to be close, but does NOT
        # substitute for a GCI/Richardson extrapolation and statistical CI.
        drift = abs(values[-1] - values[-2]) / max(abs(values[-1]), 1.0e-12)
        check(drift <= MAX_MESH_QOI_DRIFT,
              f"{scenario}: {name} still mesh-sensitive ({drift:.3%} > 5%)")
        stability[name] = drift
        # Compute a three-grid observed order / fine-grid GCI only when the
        # requested refinement ratio is sufficiently uniform and differences
        # are monotone. Do not invent p from non-asymptotic or noisy values.
        h1, h2, h3 = (m["characteristic_h"] for m in meshes[-3:])
        r12, r23 = h1 / h2, h2 / h3
        d12 = values[-3] - values[-2]
        d23 = values[-2] - values[-1]
        if abs(d12) <= 1e-12 and abs(d23) <= 1e-12:
            observed[name] = {"status": "ROUND_OFF_PLATEAU_GCI_UNDEFINED"}
        elif (abs(r12 / r23 - 1.0) <= 0.05 and d12 * d23 > 0.0
              and abs(d12) > abs(d23)):
            order = math.log(abs(d12 / d23)) / math.log(r23)
            check(0.5 <= order <= 5.0,
                  f"{scenario}: {name} observed spatial order out of bounds: {order:.3g}")
            gci = (1.25 * abs(d23) /
                   (max(abs(values[-1]), 1e-12) * (r23 ** order - 1.0)))
            check(gci <= MAX_MESH_QOI_DRIFT,
                  f"{scenario}: {name} three-grid GCI={gci:.3%} > 5%")
            observed[name] = {"status": "GCI_ESTIMATED", "order": order,
                              "gci_fine_relative": gci}
        else:
            raise EvidenceError(
                f"{scenario}: {name} three-grid GCI not demonstrable "
                "(nonuniform refinement, oscillating trend, or no asymptotic order)"
            )
    return {"scenario": scenario, "meshes": meshes,
            "finest_grid_qoi": {key: vals[-1] for key, vals in measurements.items()},
            "finest_grid_relative_mesh_drift": stability,
            "profile_mesh_convergence": profile_stability,
            "gci_evidence": observed}


def qualify_one(case: str, evidence_dir: Path, *,
                require_native_evidence: bool = False) -> dict[str, Any]:
    entry = evidence_dir / f"{case}.json"
    try:
        check(entry.is_file(), f"{case}: no CFDX solver evidence file ({entry})")
        data = json.loads(entry.read_text(encoding="utf-8"))
        check(data.get("schema_version") == 1 and data.get("case_id") == case,
              f"{case}: invalid schema_version/case_id")
        check(data.get("solver") == "CFDX" and data.get("production_run") is True,
              f"{case}: only real CFDX production-run outputs accepted")
        check(isinstance(data.get("solver_revision"), str) and
              re.fullmatch(r"[0-9a-f]{40}", data["solver_revision"]) is not None,
              f"{case}: exact 40-character git revision required")
        model = data.get("model")
        check(model in ("RANS_SST", "RSM", "LES", "DNS"),
              f"{case}: unsupported or unrecorded turbulence closure")
        check(isinstance(data.get("case_configuration_sha256"), str) and
              re.fullmatch(r"[0-9a-f]{64}", data["case_configuration_sha256"]) is not None,
              f"{case}: exact case configuration hash required")
        checked_file(evidence_dir, data.get("case_configuration_file", ""),
                     data["case_configuration_sha256"])
        checked_file(evidence_dir, data.get("solver_log_file", ""),
                     data.get("solver_log_sha256", ""))
        scenarios = data.get("scenarios")
        check(isinstance(scenarios, list), f"{case}: scenarios missing")
        found = {s.get("id") for s in scenarios if isinstance(s, dict)}
        required = set(SCENARIOS[case])
        check(found == required and len(scenarios) == len(required),
              f"{case}: complete matrix required, missing={sorted(required - found)}, "
              f"extra={sorted(found - required)}")
        if case in ("TH-01", "TH-07", "TH-08") and model in ("LES", "DNS"):
            # Resolved turbulence needs sampling/stationarity gates; until
            # dedicated statistically independent block records exist, we
            # do not let the mean-field checks certify DNS/LES quality.
            check(data.get("statistics_converged") is True and
                  type(data.get("independent_averaging_blocks")) is int and
                  data["independent_averaging_blocks"] >= 2,
                  f"{case}: converged independent time-averaging blocks required")
        evaluated = [check_meshes(case, scenario["id"], scenario.get("runs"), evidence_dir,
                                 require_native_evidence=require_native_evidence)
                     for scenario in scenarios]
        return {"case_id": case, "status": "PREQUALIFICATION_PASS",
                "physical_validation": "NOT_FULLY_QUALIFIED_MODEL_SCOPE_AND_GCI_PENDING",
                "model": model, "solver_revision": data["solver_revision"],
                "scenarios": evaluated}
    except (EvidenceError, OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        return {"case_id": case, "status": "INCOMPLETE_OR_FAIL",
                "physical_validation": "NOT_QUALIFIED", "reason": str(exc)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-pass", action="store_true",
                        help="Fail unless all cases pass, with native CFDX trace required")
    parser.add_argument("--require-native-evidence", action="store_true",
                        help="Require hash-verified CFDX run_evidence.json and native outputs")
    args = parser.parse_args()
    cases = list(dict.fromkeys(args.case or CASES))
    # The blocking qualification path must not accept manually asserted
    # solver_executed/converged flags without authentic solver output artifacts.
    native = args.require_native_evidence or args.require_pass
    results = [qualify_one(case, args.evidence_dir.resolve(),
                           require_native_evidence=native) for case in cases]
    complete = all(r["status"] == "PREQUALIFICATION_PASS" for r in results)
    report = {"schema_version": 1, "scope": "independent CFDX output prequalification",
              "never_a_full_physical_validation": True,
              "native_solver_artifacts_required": native,
              "status": "PREQUALIFICATION_ONLY" if complete else "INCOMPLETE_OR_FAIL",
              "cases": results}
    output = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8")
    print(output, end="")
    return 2 if args.require_pass and not complete else 0


if __name__ == "__main__":
    raise SystemExit(main())
