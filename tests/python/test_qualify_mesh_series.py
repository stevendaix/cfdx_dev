"""Negative and positive self-tests, never physical CFD qualification evidence."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts" / "qualify_mesh_series.py"
SHA = "a" * 40  # Synthetic fixture revision, NOT a CFDX validation claim.


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_fixture(root: Path, physics: str = "naca_laminar") -> dict:
    naca = physics.startswith("naca")
    names = ("cd_pressure", "cd_viscous", "cd_total", "cl") if naca else (
        ("cf", "delta99") if physics == "sst" else ("heat_flux", "temperature_mean")
    )
    if naca:
        reference = {"cd_pressure": (0.02, 0.02, 0.02),
                     "cd_viscous": (0.04, 0.04, 0.02),
                     "cd_total": (0.06, 0.06, 0.02),
                     "cl": (0.0, 0.001, 0.02)}
        amplitude = {"cd_pressure": 0.02, "cd_viscous": 0.01,
                     "cd_total": 0.03, "cl": 0.002}
    else:
        reference = {name: (float(i + 1), float(i + 1), 0.02)
                     for i, name in enumerate(names)}
        amplitude = {name: 0.02 for name in names}
    criteria = {
        "max_true_residual": 1e-8,
        "min_refinement_ratio": 1.5,
        "max_refinement_ratio_difference": 0.01,
        "min_observed_order": 1.5,
        "max_gci": 0.02,
    }
    balance_keys = ("energy_global", "energy_local") if physics == "thermal" else (
        "mass_global", "mass_local"
    )
    criteria.update({"max_" + key: 1e-7 for key in balance_keys})
    if naca:
        criteria["max_drag_decomposition_error"] = 1e-12
    case = {
        "schema_version": "1.0",
        "case_id": "synthetic-only",
        "physics": physics,
        "source_sha": SHA,
        "criteria": criteria,
        "references": {name: {"value": v[0], "scale": v[1],
                              "max_relative_error": v[2]}
                       for name, v in reference.items()},
        "levels": [],
    }
    for i, h in enumerate((0.25, 0.125, 0.0625)):
        mesh = root / f"mesh{i}.msh"
        mesh.write_text(f"synthetic mesh {i} with characteristic h={h}\n", encoding="utf-8")
        qoi = {name: reference[name][0] + amplitude[name] * h * h
               for name in names}
        result = {
            "case_id": "synthetic-only",
            "source_sha": SHA,
            "physics": physics,
            "status": "COMPLETED",
            "converged": True,
            "mesh_sha256": sha256(mesh),
            "cells": (16, 64, 256)[i],
            "h": h,
            "true_residual": 1e-10,
            "balances": {key: 1e-10 for key in balance_keys},
            "qoi": qoi,
        }
        rp = root / f"run{i}.json"
        rp.write_text(json.dumps(result, sort_keys=True), encoding="utf-8")
        case["levels"].append({
            "mesh": {"path": mesh.name, "sha256": sha256(mesh)},
            "result": {"path": rp.name, "sha256": sha256(rp)},
        })
    return case


def call(root: Path, case: dict, sha: str = SHA):
    path = root / "case.json"
    path.write_text(json.dumps(case, sort_keys=True), encoding="utf-8")
    report = root / "verdict.json"
    proc = subprocess.run([sys.executable, str(RUNNER), "--input", str(path),
                           "--expected-source-sha", sha, "--report", str(report)],
                          capture_output=True, text=True, check=False)
    return proc, json.loads(report.read_text(encoding="utf-8"))


class ThreeMeshEvidenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.case = make_fixture(self.root)

    def fail_with(self, expected: str) -> None:
        result, report = call(self.root, self.case)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(report["status"], "FAIL")
        self.assertIn(expected, report["reason"])

    def test_positive_synthetic_naca(self) -> None:
        result, report = call(self.root, self.case)
        self.assertEqual(result.returncode, 0, report)
        self.assertEqual(report["status"], "PASS")
        self.assertAlmostEqual(report["metrics"]["cd_total"]["observed_order"], 2.0)
        self.assertEqual(len(report["artifacts"]), 3)

    def test_positive_synthetic_sst_and_thermal(self) -> None:
        for physics in ("sst", "thermal", "naca_sst", "laminar"):
            case = make_fixture(self.root, physics)
            result, report = call(self.root, case)
            self.assertEqual(result.returncode, 0, (physics, report))

    def test_reject_missing_grid(self) -> None:
        self.case["levels"].pop()
        self.fail_with("exactly three")

    def test_reject_wrong_revision(self) -> None:
        result, report = call(self.root, self.case, "b" * 40)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("source SHA", report["reason"])

    def test_reject_tampered_mesh(self) -> None:
        (self.root / "mesh1.msh").write_text("tampered after solver run")
        self.fail_with("SHA-256 mismatch")

    def test_reject_tampered_result(self) -> None:
        data = json.loads((self.root / "run1.json").read_text())
        data["converged"] = False
        (self.root / "run1.json").write_text(json.dumps(data))
        self.fail_with("SHA-256 mismatch")

    def test_reject_solver_nonconvergence(self) -> None:
        self.edit_run(1, converged=False)
        self.fail_with("did not converge")

    def test_reject_incomplete_mass_conservation(self) -> None:
        self.edit_run(2, balances={"mass_global": 1e-10})
        self.fail_with("mass_local")

    def test_reject_large_true_residual(self) -> None:
        self.edit_run(0, true_residual=0.1)
        self.fail_with("true linear residual")

    def test_reject_nonmonotone_physical_series(self) -> None:
        data = json.loads((self.root / "run1.json").read_text())
        data["qoi"]["cd_total"] = 0.059
        self.edit_run(1, qoi=data["qoi"])
        self.fail_with("drag does not sum")
        # An independently conserved but divergent drag series must also fail.
        data["qoi"]["cd_pressure"] = 0.019
        data["qoi"]["cd_viscous"] = 0.04
        self.edit_run(1, qoi=data["qoi"])
        self.fail_with("non-monotone")

    def test_reject_fake_total_drag(self) -> None:
        data = json.loads((self.root / "run2.json").read_text())
        data["qoi"]["cd_total"] += 0.01
        self.edit_run(2, qoi=data["qoi"])
        self.fail_with("drag does not sum")

    def test_reject_nonfinite_json(self) -> None:
        path = self.root / "run2.json"
        text = path.read_text()
        path.write_text(text.replace('"true_residual": 1e-10', '"true_residual": NaN'))
        self.case["levels"][2]["result"]["sha256"] = sha256(path)
        self.fail_with("non-finite JSON")

    def test_reject_missing_reference_component(self) -> None:
        del self.case["references"]["cd_viscous"]
        self.fail_with("NACA needs separate")

    def test_reject_escaped_artifact(self) -> None:
        self.case["levels"][0]["mesh"]["path"] = "../not-the-mesh.msh"
        self.fail_with("missing, empty or out-of-tree")

    def edit_run(self, index: int, **fields) -> None:
        path = self.root / f"run{index}.json"
        data = json.loads(path.read_text())
        data.update(fields)
        path.write_text(json.dumps(data, sort_keys=True))
        self.case["levels"][index]["result"]["sha256"] = sha256(path)


if __name__ == "__main__":
    unittest.main()
