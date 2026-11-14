"""V78 geometry-comparability checks: dimensions are necessary, not sufficient."""
import importlib
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
gate_reference_geometry = importlib.import_module(
    "scripts.bfs_geometry_comparability_v78"
).gate_reference_geometry


class TestBfsGeometryComparabilityV78(unittest.TestCase):
    def test_reference_six_block_domain(self):
        result = gate_reference_geometry({
            "h_m": 0.0127, "h_upstream": 8.0, "h_downstream": 9.0,
            "x_min_h": -130.0, "x_max_h": 50.0,
        })
        self.assertTrue(result["pass_reference_domain"])
        self.assertFalse(result["physical_validation_authorized"])

    def test_shortened_pilot_domain_rejected(self):
        result = gate_reference_geometry({
            "h_m": 0.0127, "h_upstream": 7.0, "h_downstream": 8.0,
            "x_min_h": -6.0, "x_max_h": 40.0,
        })
        self.assertFalse(result["pass_reference_domain"])

    def test_small_height_error_rejected(self):
        result = gate_reference_geometry({
            "h_m": 0.0127, "h_upstream": 7.999, "h_downstream": 9.0,
            "x_min_h": -130.0, "x_max_h": 50.0,
        })
        self.assertFalse(result["pass_reference_domain"])

    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            gate_reference_geometry({
                "h_m": 0.0127, "h_upstream": float("nan"),
                "h_downstream": 9.0, "x_min_h": -130.0, "x_max_h": 50.0,
            })


if __name__ == "__main__":
    unittest.main()
