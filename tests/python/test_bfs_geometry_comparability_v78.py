"""V78 geometry-comparability checks: dimensions are necessary, not sufficient."""

import importlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
v78 = importlib.import_module("scripts.bfs_geometry_comparability_v78")
gate_reference_geometry = v78.gate_reference_geometry
pressure_failure_evidence = v78.pressure_failure_evidence


class TestBfsGeometryComparabilityV78(unittest.TestCase):
    def test_reference_six_block_domain(self):
        result = gate_reference_geometry(
            {
                "h_m": 0.0127,
                "h_upstream": 8.0,
                "h_downstream": 9.0,
                "x_min_h": -130.0,
                "x_max_h": 50.0,
            }
        )
        self.assertTrue(result["pass_reference_domain"])
        self.assertFalse(result["physical_validation_authorized"])

    def test_shortened_pilot_domain_rejected(self):
        result = gate_reference_geometry(
            {
                "h_m": 0.0127,
                "h_upstream": 7.0,
                "h_downstream": 8.0,
                "x_min_h": -6.0,
                "x_max_h": 40.0,
            }
        )
        self.assertFalse(result["pass_reference_domain"])

    def test_small_height_error_rejected(self):
        result = gate_reference_geometry(
            {
                "h_m": 0.0127,
                "h_upstream": 7.999,
                "h_downstream": 9.0,
                "x_min_h": -130.0,
                "x_max_h": 50.0,
            }
        )
        self.assertFalse(result["pass_reference_domain"])

    def test_native_pressure_rejections_are_extracted_not_summary_duplicated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            log = root / "native.log"
            execution = root / "execution.json"
            rejected = (
                "CFDX nonlinear iteration 1 rejected; retry=1 reason=solve_steady_incompressible: "
                "pressure-correction solve did not converge "
                "(status=1, iterations=1000, residual=1e-12, relative=5e-10)"
            )
            log.write_text(
                rejected
                + "\n"
                + rejected.replace("retry=1", "retry=2")
                + "\nCFDX production solver error: pressure-correction solve did not converge "
                "(status=1, iterations=1000, residual=1e-12, relative=5e-10)\n"
            )
            execution.write_text(
                json.dumps({"converged": False, "process_exit_code": 2})
            )
            evidence = pressure_failure_evidence(log, execution)
            self.assertEqual(len(evidence["trials"]), 2)
            self.assertFalse(evidence["physical_validation_authorized"])
            execution.write_text(
                json.dumps({"converged": True, "process_exit_code": 0})
            )
            with self.assertRaises(ValueError):
                pressure_failure_evidence(log, execution)

    def test_pressure_summary_without_native_rejection_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            log = root / "native.log"
            execution = root / "execution.json"
            log.write_text(
                "CFDX production solver error: pressure-correction solve did not converge "
                "(status=1, iterations=1000, residual=1e-12, relative=5e-10)\n"
            )
            execution.write_text(
                json.dumps({"converged": False, "process_exit_code": 2})
            )
            with self.assertRaisesRegex(ValueError, "no native pressure-rejection"):
                pressure_failure_evidence(log, execution)

    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            gate_reference_geometry(
                {
                    "h_m": 0.0127,
                    "h_upstream": float("nan"),
                    "h_downstream": 9.0,
                    "x_min_h": -130.0,
                    "x_max_h": 50.0,
                }
            )


if __name__ == "__main__":
    unittest.main()
