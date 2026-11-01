from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "validation"
    / "analyze_naca_pressure_spectrum.py"
)
SPEC = importlib.util.spec_from_file_location("analyze_naca_pressure_spectrum", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class NacaPressureSpectrumTests(unittest.TestCase):
    def test_condition_estimate_requires_positive_finite_eigenvalues(self):
        self.assertEqual(MODULE.condition_estimate(2.0, 10.0), 5.0)
        self.assertIsNone(MODULE.condition_estimate(0.0, 10.0))
        self.assertIsNone(MODULE.condition_estimate(-1.0, 10.0))
        self.assertIsNone(MODULE.condition_estimate(float("nan"), 10.0))

    def test_parse_rhs_csv_preserves_physical_cell_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.rhs.csv"
            path.write_text(
                "row,cell,rhs\n0,1,-2.5\n1,3,4.0\n",
                encoding="utf-8",
            )
            records = MODULE.parse_rhs_csv(path)
        self.assertEqual(records[0]["cell"], 1)
        self.assertEqual(records[1]["row"], 1)
        self.assertEqual(records[1]["rhs"], 4.0)


if __name__ == "__main__":
    unittest.main()
