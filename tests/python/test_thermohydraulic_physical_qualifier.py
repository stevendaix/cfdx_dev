"""Mutation/coverage tests: synthetic fixtures exercise the checker, NOT CFDX."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import qualify_thermohydraulic_benchmarks as qual  # noqa: E402


def pinned(folder: Path, name: str, body: str) -> tuple[str, str]:
    (folder / name).write_text(body, encoding="utf-8")
    return name, hashlib.sha256(body.encode("utf-8")).hexdigest()


def make_th03(folder: Path) -> dict:
    config, config_hash = pinned(folder, "input.toml", "# synthetic, not a CFDX case\n")
    log, log_hash = pinned(folder, "stdout.log", "CFDX synthetic fixture\n")
    reference = "x,value,sigma\n" + "".join(f"{i},{1+i/10:.8f},0.01\n" for i in range(6))
    simulation = "x,value\n" + "".join(f"{i},{1+i/10:.8f}\n" for i in range(6))
    ref_path, ref_hash = pinned(folder, "reference.csv", reference)
    sim_path, sim_hash = pinned(folder, "simulation.csv", simulation)
    profiles = []
    for name, station in (("cf", "lower_wall"), ("cp", "lower_wall"),
                          *(("u_mean", s) for s in
                            ("x/h=-4", "x/h=1", "x/h=4", "x/h=6", "x/h=10"))):
        profiles.append({"name": name, "station": station,
                         "units": qual.PROFILE_UNITS["TH-03"][name],
                         "reference_source_url": "https://example.org/synthetic-fixture-not-real",
                         "reference_file": ref_path, "reference_sha256": ref_hash,
                         "simulation_file": sim_path, "simulation_sha256": sim_hash})
    runs = []
    for i, (qoi, cells, h) in enumerate(((6.00, 200, 0.25),
                                        (6.20, 800, 0.125),
                                        (6.26, 3200, 0.0625))):
        runs.append({"mesh_id": f"level{i}", "cells": cells, "characteristic_h": h,
                     "solver_executed": True, "converged": True,
                     "true_residual_max": 1e-10, "mass_imbalance_relative": 1e-7,
                     "qoi": {"reattachment_x_over_h": qoi},
                     "profiles": copy.deepcopy(profiles)})
    data = {"schema_version": 1, "case_id": "TH-03", "solver": "CFDX",
            "model": "RANS_SST", "production_run": True,
            "solver_revision": "a" * 40, "case_configuration_file": config,
            "case_configuration_sha256": config_hash,
            "solver_log_file": log, "solver_log_sha256": log_hash,
            "scenarios": [{"id": "Re_h_37500", "runs": runs}]}
    (folder / "TH-03.json").write_text(json.dumps(data), encoding="utf-8")
    return data


def write_manifest(folder: Path, data: dict) -> None:
    (folder / "TH-03.json").write_text(json.dumps(data), encoding="utf-8")


def test_missing_all_four_fails_closed(tmp_path: Path) -> None:
    report = [qual.qualify_one(case, tmp_path) for case in qual.CASES]
    assert len(report) == 4
    assert all(r["status"] == "INCOMPLETE_OR_FAIL" for r in report)
    assert all(r["physical_validation"] == "NOT_QUALIFIED" for r in report)


def test_synthetic_fixture_exercises_prequalification_not_full_validation(tmp_path: Path) -> None:
    make_th03(tmp_path)
    result = qual.qualify_one("TH-03", tmp_path)
    assert result["status"] == "PREQUALIFICATION_PASS", result
    assert result["physical_validation"] == "NOT_FULLY_QUALIFIED_MODEL_SCOPE_AND_GCI_PENDING"
    assert result["scenarios"][0]["gci_evidence"]["reattachment_x_over_h"]["order"] > 1


@pytest.mark.parametrize("mutate,expected", [
    (lambda d: d.update({"production_run": False}), "production"),
    (lambda d: d["scenarios"][0]["runs"].pop(), "three distinct"),
    (lambda d: d["scenarios"][0]["runs"][-1]["qoi"].update(
        {"reattachment_x_over_h": 7.2}), "outside reference"),
    (lambda d: d["scenarios"][0]["runs"][-1].update(
        {"mass_imbalance_relative": 0.2}), "mass imbalance"),
    (lambda d: d["scenarios"][0]["runs"][-1].update(
        {"true_residual_max": float("nan")}), "nonfinite"),
    (lambda d: d["scenarios"][0]["runs"][0].update({"converged": False}), "converged"),
    (lambda d: d["scenarios"][0]["runs"][-1].update({"characteristic_h": 0.12}), "refinement"),
    (lambda d: d["scenarios"][0]["runs"][-1]["profiles"][0].update(
        {"reference_sha256": "0" * 64}), "SHA-256 mismatch"),
    (lambda d: d["scenarios"][0]["runs"][-1]["profiles"][0].update(
        {"units": "m/s"}), "normalization"),
    (lambda d: d["scenarios"][0]["runs"][-1]["profiles"].pop(), "missing BFS velocity stations"),
    (lambda d: d.update({"solver_log_sha256": "1" * 64}), "SHA-256 mismatch"),
])
def test_mutations_must_fail(tmp_path: Path, mutate, expected: str) -> None:
    data = make_th03(tmp_path)
    mutate(data)
    write_manifest(tmp_path, data)
    result = qual.qualify_one("TH-03", tmp_path)
    assert result["status"] == "INCOMPLETE_OR_FAIL", result
    assert expected in result["reason"], result


def test_profile_interpolation_rejects_missing_domain(tmp_path: Path) -> None:
    data = make_th03(tmp_path)
    p = data["scenarios"][0]["runs"][0]["profiles"][0]
    changed = "x,value\n" + "".join(f"{i},{1+i/10:.8f}\n" for i in range(1, 7))
    _, checksum = pinned(tmp_path, "simulation.csv", changed)
    for run in data["scenarios"][0]["runs"]:
        for profile in run["profiles"]:
            profile["simulation_sha256"] = checksum
    write_manifest(tmp_path, data)
    result = qual.qualify_one("TH-03", tmp_path)
    assert "no extrapolation" in result["reason"]


def test_all_scenarios_required_and_unique(tmp_path: Path) -> None:
    data = make_th03(tmp_path)
    data["scenarios"].append(copy.deepcopy(data["scenarios"][0]))
    write_manifest(tmp_path, data)
    assert "complete matrix required" in qual.qualify_one("TH-03", tmp_path)["reason"]


def test_profiles_must_improve_on_three_meshes_not_just_scalar_qoi(tmp_path: Path) -> None:
    data = make_th03(tmp_path)
    # All profiles still lie inside independent reference uncertainty, but their
    # mesh-to-mesh difference grows. An integral-only GCI gate would miss this.
    changed = "x,value\n" + "".join(f"{i},{1+i/10+0.002:.8f}\n" for i in range(6))
    _, digest = pinned(tmp_path, "fine-only.csv", changed)
    for profile in data["scenarios"][0]["runs"][-1]["profiles"]:
        profile["simulation_file"] = "fine-only.csv"
        profile["simulation_sha256"] = digest
    write_manifest(tmp_path, data)
    result = qual.qualify_one("TH-03", tmp_path)
    assert result["status"] == "INCOMPLETE_OR_FAIL", result
    assert "does not improve with refinement" in result["reason"]


def test_reference_dataset_must_not_change_between_meshes(tmp_path: Path) -> None:
    data = make_th03(tmp_path)
    changed = "x,value,sigma\n" + "".join(f"{i},{1+i/10:.8f},0.02\n" for i in range(6))
    _, digest = pinned(tmp_path, "other-reference.csv", changed)
    for profile in data["scenarios"][0]["runs"][-1]["profiles"]:
        profile["reference_file"] = "other-reference.csv"
        profile["reference_sha256"] = digest
    write_manifest(tmp_path, data)
    result = qual.qualify_one("TH-03", tmp_path)
    assert result["status"] == "INCOMPLETE_OR_FAIL", result
    assert "reference dataset" in result["reason"]


def provide_synthetic_native_traces(folder: Path, data: dict) -> None:
    """Test fixtures only: exercise the strict artifact-consistency guard."""
    for run in data["scenarios"][0]["runs"]:
        native_dir = folder / "native" / run["mesh_id"]
        native_dir.mkdir(parents=True)
        pinned(native_dir, "solver.log", "fixture log\n")
        files = {
            "execution.json": {"format":"CFDX-EXECUTION", "process_exit_code":0,
                               "converged":True},
            "convergence.json": {"format":"CFDX-CONVERGENCE", "converged":True,
                                 "history":[{"momentum_equation_residual_rhs_relative":1e-10,
                                             "continuity_normalized":1e-11}]},
            "benchmark_observables/mass_conservation.json": {
                "format":"CFDX-THERMOHYDRAULIC-CONSERVATION",
                "physical_qualification":"NOT_QUALIFIED", "cells":run["cells"],
                "global_mass_imbalance_relative":1e-7,
                "maximum_cell_mass_imbalance_relative":1e-8},
        }
        checksums = {}
        for file, content in files.items():
            path = native_dir / file
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(json.dumps(content), encoding="utf-8")
            checksums[file] = hashlib.sha256(path.read_bytes()).hexdigest()
        for file in ("benchmark_observables/cell_fields.csv",
                     "benchmark_observables/face_fluxes.csv"):
            path = native_dir/file
            path.write_text("fixture,not_real\n", encoding="utf-8")
            checksums[file] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest = {"case_id":"TH-03", "scenario_id":"Re_h_37500",
                    "mesh_id":run["mesh_id"], "status":"EXECUTED_UNQUALIFIED",
                    "physical_qualification":"NOT_QUALIFIED", "solver_executed":True,
                    "exit_code":0, "solver_binary_sha256":"a"*64,
                    "case_file_sha256":"b"*64,
                    "solver_log_sha256":hashlib.sha256((native_dir/"solver.log").read_bytes()).hexdigest(),
                    "true_residual_max":1e-10,"mass_imbalance_relative":1e-7,
                    "files_sha256":checksums}
        trace_path = native_dir / "run_evidence.json"
        trace_path.write_text(json.dumps(manifest), encoding="utf-8")
        run["native_run_evidence_file"] = str(trace_path.relative_to(folder))
        run["native_run_evidence_sha256"] = hashlib.sha256(trace_path.read_bytes()).hexdigest()
    write_manifest(folder, data)


def test_strict_native_evidence_required_even_for_synthetic_prequalification(tmp_path: Path) -> None:
    data = make_th03(tmp_path)
    assert qual.qualify_one("TH-03",tmp_path)["status"] == "PREQUALIFICATION_PASS"
    result = qual.qualify_one("TH-03",tmp_path,require_native_evidence=True)
    assert result["status"] == "INCOMPLETE_OR_FAIL"
    assert "file path missing" in result["reason"]
    provide_synthetic_native_traces(tmp_path,data)
    # Merely a contract unit test, NOT an actual CFDX run or physical evidence.
    assert qual.qualify_one("TH-03",tmp_path,require_native_evidence=True)["status"] == "PREQUALIFICATION_PASS"
    data["scenarios"][0]["runs"][0]["mass_imbalance_relative"] = 3e-7
    write_manifest(tmp_path, data)
    result = qual.qualify_one("TH-03",tmp_path,require_native_evidence=True)
    assert result["status"] == "INCOMPLETE_OR_FAIL"
    assert "not backed by the native CFDX outputs" in result["reason"]


def test_strict_native_evidence_detects_modified_solver_output(tmp_path: Path) -> None:
    data = make_th03(tmp_path)
    provide_synthetic_native_traces(tmp_path,data)
    (tmp_path / "native" / "level0" / "convergence.json").write_text(
        '{"format":"CFDX-CONVERGENCE","converged":false}',encoding="utf-8")
    result = qual.qualify_one("TH-03",tmp_path,require_native_evidence=True)
    assert result["status"] == "INCOMPLETE_OR_FAIL"
    assert "SHA-256 mismatch" in result["reason"]
