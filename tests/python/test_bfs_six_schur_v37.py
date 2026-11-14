"""V37 six-Schur physical-input equality and fail-closed evidence checks."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from run_bfs_six_schur_v37 import (  # noqa: E402
    PRECONDITIONER, SCHUR_MODELS, physical_input_identity, prepare_cases,
    run_one, summarize_pcd_nested_fp,
)


class TestBFSSixSchurV37(unittest.TestCase):
    def test_six_identical_physical_inputs_with_six_explicit_schur(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = prepare_cases(Path(tmp) / "work", iterations=1,
                                   reference_case=None, upstream=4,
                                   downstream=12, ny=10, step_cells=2,
                                   wall_stretch=0.0)
            self.assertEqual(set(result["cases"]), set(SCHUR_MODELS))
            self.assertEqual(len({entry["physical_sha256"]
                                  for entry in result["cases"].values()}), 1)
            self.assertEqual(len({entry["case_sha256"]
                                  for entry in result["cases"].values()}), 6)
            self.assertTrue(all(entry["n_cells"] == 152
                                for entry in result["cases"].values()))
            self.assertEqual(result["preconditioner_policy"], PRECONDITIONER)
            self.assertNotIn("PASS", result["physical_validation"])
            # Every case retains the same turbulence, boundary and initial data.
            for schur, entry in result["cases"].items():
                with h5py.File(entry["case"], "r") as h5:
                    setup = json.loads(h5.attrs["case_setup_json"])
                    selection = {x["family"]: x["configuration_key"]
                                 for x in setup["numerics"]["selection"]["entries"]}
                    self.assertEqual(selection["schur"], f"schur.{schur}")
                    self.assertEqual(selection["preconditioner"],
                                     "preconditioner.coupled_block_schur")
                    self.assertEqual(set(setup["initial_condition"]["scalar_fields"]),
                                     {"k", "omega"})

    def test_forensic_operator_capture_enabled_identically_for_all_six(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = prepare_cases(Path(tmp) / "work", iterations=1,
                                   reference_case=None, upstream=4,
                                   downstream=12, ny=10, step_cells=2,
                                   wall_stretch=0.0, capture_operator=True)
            self.assertTrue(result["operator_capture_requested"])
            self.assertEqual(len({entry["physical_sha256"] for entry
                                  in result["cases"].values()}), 1)
            for entry in result["cases"].values():
                with h5py.File(entry["case"], "r") as h5:
                    setup = json.loads(h5.attrs["case_setup_json"])
                    self.assertEqual(
                        setup["numerics"]["raw_settings"]["coupled_diagnostic_capture_prefix_v16"],
                        result["operator_capture_prefix"])

    def test_changed_input_is_detected_even_with_unchanged_mesh_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = prepare_cases(Path(tmp) / "work", iterations=2,
                                   reference_case=None, upstream=4,
                                   downstream=12, ny=10, step_cells=2,
                                   wall_stretch=0.0)
            case = Path(result["cases"]["lsc"]["case"])
            with h5py.File(case, "r+") as h5:
                setup = json.loads(h5.attrs["case_setup_json"])
                setup["boundary_conditions"][0]["velocity_vector"][0] = 39.0
                h5.attrs["case_setup_json"] = json.dumps(setup)
            modified = physical_input_identity(case)
            self.assertNotEqual(modified["physical_sha256"],
                                result["common_physical_sha256"])
            self.assertEqual(modified["topology_hash"],
                             result["cases"]["lsc"]["topology_hash"])

    def test_preconditioner_replacement_cannot_pass_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = prepare_cases(Path(tmp) / "work", iterations=2,
                                   reference_case=None, upstream=4,
                                   downstream=12, ny=10, step_cells=2,
                                   wall_stretch=0.0)
            case = Path(result["cases"]["pcd"]["case"])
            with h5py.File(case, "r+") as h5:
                setup = json.loads(h5.attrs["case_setup_json"])
                for entry in setup["numerics"]["selection"]["entries"]:
                    if entry["family"] == "preconditioner":
                        entry["configuration_key"] = "preconditioner.pcd"
                h5.attrs["case_setup_json"] = json.dumps(setup)
            with self.assertRaisesRegex(ValueError, "fixed solver policy"):
                physical_input_identity(case)

    def test_fails_on_existing_output_to_prevent_stale_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "work"
            out.mkdir()
            (out / "old_run.json").write_text('{"converged":true}')
            with self.assertRaises(FileExistsError):
                prepare_cases(out, iterations=1, reference_case=None,
                              upstream=4, downstream=12, ny=10, step_cells=2,
                              wall_stretch=0.0)

    def test_nested_pcd_diagnostic_separates_quality_from_acceptance(self):
        record = """COUPLED_NESTED_PCD_FP status=1 iterations=100 original_true_relative_rhs=0.006 accepted=1
COUPLED_NESTED_PCD_FP status=0 iterations=32 original_true_relative_rhs=1e-11 accepted=1
COUPLED_NESTED_PCD_FP status=1 iterations=100 original_true_relative_rhs=nan accepted=1
"""
        result = summarize_pcd_nested_fp(record)
        self.assertEqual(result["calls"], 2)
        self.assertEqual(result["malformed_or_nonfinite"], 1)
        self.assertEqual(result["status_max_iter"], 1)
        self.assertEqual(result["accepted_actions"], 2)
        self.assertEqual(result["original_true_relative_over_1e3"], 1)
        self.assertEqual(result["original_true_relative_under_1e8"], 1)
        self.assertFalse(result["physical_qualified"])

    def test_nonexistent_solver_never_makes_a_false_certificate(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "work"
            result = prepare_cases(out, iterations=1, reference_case=None,
                                   upstream=4, downstream=12, ny=10,
                                   step_cells=2, wall_stretch=0.0)
            case = {**result["cases"]["pcd"], "iterations": 1}
            outcome = run_one(Path(tmp) / "missing_executable", out,
                              "pcd", case, wall_timeout=1)
            self.assertEqual(outcome["execution"], "FAILED_TO_START")
            self.assertIsNone(outcome["exit_code"])
            self.assertFalse(outcome["nonlinear_convergence_proven"])
            self.assertNotEqual(outcome["v8_evidence"]["physical_validation"], "PASS")


if __name__ == "__main__":
    unittest.main()
