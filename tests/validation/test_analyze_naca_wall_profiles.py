"""Tests for independent symmetry/integral audit of NACA wall profiles."""

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


PATH = Path(__file__).resolve().parents[2] / "scripts/analyze_naca_wall_profiles.py"
SPEC = importlib.util.spec_from_file_location("wall_profile_audit", PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def fake_profiles(path):
    # 16 exact mirrored pairs, no fabricated NACA aerodynamic benchmark.
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(MODULE.FIELDS)
        for i in range(16):
            for side in (-1, 1):
                writer.writerow((2*i+(side+1)//2, i/15, side*0.04,
                                 -0.1+i*0.01, 0.03+i*0.001,
                                 0.005, 0.003))


class ProfileAuditTests(unittest.TestCase):
    def test_symmetric_profiles_and_integrals(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wall.csv"
            fake_profiles(path)
            report = MODULE.analyze(path, expected_pressure=0.16,
                                    expected_viscous=0.096)
            self.assertEqual(report["paired_faces_per_side"], 16)
            self.assertEqual(report["physical_reference_validation"], "NOT_PERFORMED")

    def test_bad_reported_force_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wall.csv"
            fake_profiles(path)
            with self.assertRaises(ValueError):
                MODULE.analyze(path, expected_pressure=0.17)

    def test_bad_symmetry_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wall.csv"
            fake_profiles(path)
            text = path.read_text(encoding="utf-8").replace("-0.1,0.03", "-0.9,0.03", 1)
            path.write_text(text, encoding="utf-8")
            with self.assertRaises(ValueError):
                MODULE.analyze(path)

    def test_missing_data_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wall.csv"
            path.write_text("face,Cp\n1,0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                MODULE.analyze(path)


if __name__ == "__main__":
    unittest.main()
