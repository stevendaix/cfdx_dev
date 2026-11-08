"""Mutation tests: a PASS is impossible when any independent gate is broken.

All fixture bytes are synthetic; these are checks of the verifier, not CFD runs.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import qualification_evidence_gate as gate


def file_item(path: Path, base: Path) -> dict:
    return {"path": str(path.relative_to(base)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


class EvidenceGateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.case_id = "LAM-COUETTE"
        self.case = {"id": self.case_id, "status": "PASS", "domain": "Laminaire"}
        self.dir = self.base / self.case_id
        self.dir.mkdir()
        self.setup = {"model": "laminar", "reynolds": 100.0,
                      "mach": 0.0, "alpha_degrees": 0.0}
        self.ref = {"case_id": self.case_id, "setup": dict(self.setup),
                    "citation": "Exact plane Couette u(y)=y/H",
                    "independent_of_solver": True,
                    "quantities": {
                        key: {"value": ref, "unit": "m/s" if key != "velocity_profile_l2" else "1",
                              "uncertainty": 0.0, "absolute_tolerance": 0.01,
                              "relative_tolerance": 0.0}
                        for key, ref in (("u_mean", 0.5), ("velocity_profile_l2", 0.0),
                                         ("wall_shear", 1.0))}}
        self._save(self.dir / "reference.json", self.ref)
        artifacts = {}
        for name in ("binary", "configuration", "solver_log", "oracle_source"):
            path = self.dir / f"{name}.txt"
            path.write_text(f"synthetic {name} fixture (not a solver result)\n", encoding="utf-8")
            artifacts[name] = file_item(path, self.dir)
        artifacts["reference"] = file_item(self.dir / "reference.json", self.dir)
        archive = self.dir / "source_archive.zip"
        archive.write_bytes(b"synthetic fixture only; not an actual CFDX ZIP\n")
        artifacts["source_archive"] = file_item(archive, self.dir)
        self.evidence = {"schema": gate.SCHEMA, "case_id": self.case_id,
                         "setup": self.setup,
                         "provenance": {"source_kind": "archive",
                                        "source_revision": artifacts["source_archive"]["sha256"],
                                        "artifacts": artifacts},
                         "meshes": []}
        for level, h in zip((16, 32, 64), (0.25, 0.125, 0.0625)):
            name = f"N{level}"
            residuals = {key: 1e-11 for key in
                         ("momentum", "continuity", "mass_global", "mass_local")}
            balances = {"mass": {"inward": 1.0, "outward": 1.0,
                                 "source": 0.0, "storage": 0.0,
                                 "relative_limit": 1e-8}}
            values = {"u_mean": 0.5 + 2 * h**2,
                      "velocity_profile_l2": 2 * h**2,
                      "wall_shear": 1.0 + 2 * h**2}
            rec = {"id": name, "status": "CONVERGED", "iterations": 4,
                   "runtime_seconds": 1.0, "final_scheme": "TVD",
                   "h": h, "mesh_family": "uniform-couette",
                   "mesh_quality": {"inverted_cells": 0,
                                    "minimum_cell_volume": h * h,
                                    "max_nonorthogonality_degrees": 30.0},
                   "residuals": residuals,
                   "residual_limits": {key: 1e-8 for key in residuals},
                   "independent_residuals": residuals.copy(),
                   "balances": balances, "quantities": values,
                   "artifacts": {}}
            meshpath = self.dir / f"{name}.mesh"
            meshpath.write_text(f"synthetic mesh {name}\n", encoding="utf-8")
            rec["artifacts"]["mesh"] = file_item(meshpath, self.dir)
            history = self.dir / f"{name}.csv"
            with history.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["iteration", *residuals, *values])
                writer.writeheader()
                writer.writerow({"iteration": 1, **{key: 1e-5 for key in residuals}, **values})
                for iteration in (2, 3, 4):
                    writer.writerow({"iteration": iteration, **residuals, **values})
            rec["artifacts"]["history"] = file_item(history, self.dir)
            oracle = {"schema": "cfdx-independent-fv-oracle/v1",
                      "case_id": self.case_id, "mesh_id": name,
                      "residuals": rec["independent_residuals"],
                      "balances": balances}
            oraclepath = self.dir / f"{name}.oracle.json"
            self._save(oraclepath, oracle)
            rec["artifacts"]["independent_oracle"] = file_item(oraclepath, self.dir)
            self.evidence["meshes"].append(rec)
        self.write_evidence()

    @staticmethod
    def _save(path: Path, item: dict):
        path.write_text(json.dumps(item, indent=2) + "\n", encoding="utf-8")

    def write_evidence(self):
        self._save(self.dir / "evidence.json", self.evidence)

    def result(self):
        return gate.evaluate_case(self.case, self.dir / "evidence.json")

    def assert_blocked(self, check):
        report = self.result()
        self.assertEqual(report["verification_status"], "NOT_QUALIFIED")
        self.assertIn(check, {x["check"] for x in report["checks"] if x["status"] == "FAIL"})

    def test_synthetic_positive_contract_only(self):
        report = self.result()
        self.assertEqual(report["verification_status"], "EVIDENCE_GATES_PASS")
        self.assertEqual(set(x["status"] for x in report["checks"]), {"PASS"})

    def test_registry_historical_pass_with_no_evidence_is_blocked(self):
        reg = json.loads((ROOT / "docs/validation/CFDX_QUALIFICATION_REGISTRY.json").read_text())
        report = gate.evaluate_registry(reg, self.base / "nonexistent")
        self.assertEqual(report["coverage"]["selected"], 25)
        self.assertEqual(report["coverage"]["evidence_pass"], 0)
        self.assertEqual(report["coverage"]["not_qualified"], 25)

    def test_archive_revision_must_equal_actual_archive_sha(self):
        self.evidence["provenance"]["source_revision"] = "0" * 64
        self.write_evidence()
        self.assert_blocked("provenance")

    def test_modified_source_archive_is_not_accepted(self):
        (self.dir / "source_archive.zip").write_bytes(b"altered source bytes")
        self.assert_blocked("provenance")

    def test_drifting_qoi_is_rejected_even_with_small_residuals(self):
        path = self.dir / "N32.csv"
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        rows[-3]["wall_shear"] = "1.20"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        self.evidence["meshes"][1]["artifacts"]["history"] = file_item(path, self.dir)
        self.write_evidence()
        self.assert_blocked("qoi_stability")

    def test_raw_final_qoi_cannot_be_stale(self):
        path = self.dir / "N32.csv"
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        rows[-1]["u_mean"] = "0.6"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        self.evidence["meshes"][1]["artifacts"]["history"] = file_item(path, self.dir)
        self.write_evidence()
        self.assert_blocked("qoi_stability")

    def test_domain_thermal_solid_needs_energy_not_artificial_momentum(self):
        one_mesh = {"id": "S16", "residuals": {"energy": 1e-11},
                    "residual_limits": {"energy": 1e-8},
                    "independent_residuals": {"energy": 1e-11},
                    "final_scheme": "CENTRAL"}
        report = gate.verify_residuals([one_mesh],
                    {"id": "T1-SOLID-COND", "domain": "Thermal"})
        self.assertEqual(report["max_independent_or_solver_residual"], {"energy": 1e-11})

    def test_separation_and_radiation_have_explicit_balances(self):
        flow = {"id": "S16", "balances": {}}
        with self.assertRaisesRegex(ValueError, "mass"):
            gate.verify_conservation([flow], {"id": "SEP-BFS", "domain": "Separation"})
        with self.assertRaisesRegex(ValueError, "radiation"):
            gate.verify_conservation([flow], {"id": "RAD-P1", "domain": "Radiation"})

    def test_modified_binary_is_not_accepted(self):
        (self.dir / "binary.txt").write_text("tampered", encoding="utf-8")
        self.assert_blocked("provenance")

    def test_missing_oracle_file_fails_closed(self):
        (self.dir / "N32.oracle.json").unlink()
        self.assert_blocked("raw_evidence")

    def test_modified_mesh_fails_closed(self):
        (self.dir / "N32.mesh").write_text("tampered mesh", encoding="utf-8")
        self.assert_blocked("raw_evidence")

    def test_mismatched_final_residual_in_history_fails(self):
        (self.dir / "N16.csv").write_text((self.dir / "N16.csv").read_text().replace("1e-11", "1e-3"))
        self.assert_blocked("raw_evidence")

    def test_missing_final_history_iterate_fails(self):
        path = self.dir / "N16.csv"
        path.write_text("iteration,momentum,continuity,mass_global,mass_local\n"
                        "1,1e-11,1e-11,1e-11,1e-11\n")
        rec = self.evidence["meshes"][0]
        rec["artifacts"]["history"] = file_item(path, self.dir)
        self.write_evidence()
        self.assert_blocked("raw_evidence")

    def test_physically_wrong_reference_is_not_passed(self):
        self.ref["quantities"]["u_mean"]["value"] = 0.72
        self._save(self.dir / "reference.json", self.ref)
        self.evidence["provenance"]["artifacts"]["reference"] = file_item(self.dir / "reference.json", self.dir)
        self.write_evidence()
        self.assert_blocked("reference")

    def test_reference_parameters_must_match_simulated_regime(self):
        self.ref["setup"]["reynolds"] = 200.0
        self._save(self.dir / "reference.json", self.ref)
        self.evidence["provenance"]["artifacts"]["reference"] = file_item(self.dir / "reference.json", self.dir)
        self.write_evidence()
        self.assert_blocked("case_contract")

    def test_unconverged_physics_fails_even_if_final_qoi_is_good(self):
        self.evidence["meshes"][1]["independent_residuals"]["momentum"] = 1e-3
        self.write_evidence()
        self.assert_blocked("residuals")

    def test_bad_global_mass_balance_fails(self):
        self.evidence["meshes"][2]["balances"]["mass"]["outward"] = 0.9
        self.write_evidence()
        self.assert_blocked("conservation")

    def test_wrong_grid_trend_blocks_spatial_order(self):
        self.evidence["meshes"][1]["quantities"]["u_mean"] = 0.6
        self.write_evidence()
        self.assert_blocked("spatial_order")

    def test_fine_value_matches_reference_but_gci_too_large(self):
        # Fine error=0.0078125 below tolerance=0.008, but GCI=0.009765625.
        for spec in self.ref["quantities"].values():
            spec["absolute_tolerance"] = 0.008
        self._save(self.dir / "reference.json", self.ref)
        self.evidence["provenance"]["artifacts"]["reference"] = file_item(
            self.dir / "reference.json", self.dir)
        self.write_evidence()
        self.assert_blocked("spatial_order")

    def test_identical_mesh_qoi_without_order_is_not_promoted(self):
        for mesh in self.evidence["meshes"]:
            mesh["quantities"]["u_mean"] = 0.5
            # Match all raw history observations, so only the order oracle
            # can detect this degenerate three-grid evidence.
            path = self.dir / f"{mesh['id']}.csv"
            with path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            for row in rows:
                row["u_mean"] = "0.5"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            mesh["artifacts"]["history"] = file_item(path, self.dir)
        self.write_evidence()
        self.assert_blocked("spatial_order")

    def test_unsimilar_mesh_labels_block(self):
        self.evidence["meshes"][2]["mesh_family"] = "other-grid"
        self.write_evidence()
        self.assert_blocked("mesh_series")

    def test_output_does_not_claim_physics_when_missing(self):
        output = self.base / "result.json"
        cmd = [sys.executable, str(ROOT / "scripts/qualification_evidence_gate.py"),
               "--evidence-dir", str(self.base / "none"),
               "--case", self.case_id, "--require-all", "--output", str(output)]
        p = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(p.returncode, 2)
        self.assertEqual(json.loads(output.read_text())["coverage"]["evidence_pass"], 0)


if __name__ == "__main__":
    unittest.main()
