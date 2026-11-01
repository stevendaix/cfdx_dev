from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/naca_wall_pressure_grid_audit.py"
SPEC = importlib.util.spec_from_file_location("naca_pressure_audit", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def sample(factor: int):
    return [{"face": i, "x_over_c": (i + 0.5) / (7 * factor),
             "y_over_c": 0.01, "Cp": 0.2, "Cf_tangent": 0.03,
             "Cd_pressure_face": 0.01 / factor,
             "Cd_viscous_face": 0.02 / factor}
            for i in range(7 * factor)]


class AuditTest(unittest.TestCase):
    def test_two_grid_nonqualification(self):
        # Populate every chord interval rather than implying statistical coverage.
        coarse = []
        fine = []
        for j, (lo, hi) in enumerate(MODULE.BINS):
            for k in range(1):
                coarse.append({"face": len(coarse), "x_over_c": (lo + hi) / 2,
                               "y_over_c": 0.02, "Cp": 0.4,
                               "Cf_tangent": 0.01, "Cd_pressure_face": 0.01,
                               "Cd_viscous_face": 0.02})
            for k in range(2):
                fine.append({"face": len(fine), "x_over_c": (lo + hi) / 2,
                             "y_over_c": 0.02, "Cp": 0.5,
                             "Cf_tangent": 0.01, "Cd_pressure_face": 0.005,
                             "Cd_viscous_face": 0.01})
        result = MODULE.audit(coarse, fine)
        self.assertEqual(result["gci"], "NOT_PERFORMED_NEEDS_CONVERGED_N256")
        self.assertAlmostEqual(result["cd_total_delta"], 0)

    def test_missing_region_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.audit(sample(1), sample(2))

    def test_duplicate_face_rejected(self):
        import csv
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "profile.csv"
            with p.open("w", newline="") as handle:
                w = csv.DictWriter(handle, fieldnames=MODULE.COLUMNS)
                w.writeheader()
                for row in sample(2):
                    row["face"] = 1
                    w.writerow(row)
            with self.assertRaises(ValueError):
                MODULE.read_profile(p)

    def test_nonfinite_rejected(self):
        import csv
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "profile.csv"
            with p.open("w", newline="") as handle:
                w = csv.DictWriter(handle, fieldnames=MODULE.COLUMNS)
                w.writeheader()
                for row in sample(2):
                    row["Cp"] = float("nan")
                    w.writerow(row)
            with self.assertRaises(ValueError):
                MODULE.read_profile(p)


if __name__ == "__main__":
    unittest.main()
