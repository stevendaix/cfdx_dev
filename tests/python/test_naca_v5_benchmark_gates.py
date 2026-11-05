"""Negative/positive tests of reproducible V4/V5 AMG qualification gates."""
from pathlib import Path
import sys
import unittest

BENCH = Path(__file__).resolve().parents[2] / "tools/benchmarks"
sys.path.insert(0, str(BENCH))
import compare_naca_amg_v5  # noqa: E402
import compare_naca_outer_histories  # noqa: E402


class GateTests(unittest.TestCase):
    def setUp(self):
        self.paths = [Path(f"/tmp/actual_N128_capture_{i}.manifest") for i in range(1, 4)]
        self.header = ("manifest,status,sequence,rows,iterations,true_relative_residual,"
                       "total_ms,full_setup_ms,numeric_update_ms,krylov_ms,"
                       "amg_operator_complexity,amg_levels\n")
        self.rows = [f"{p},PASS,{i},12336,9,2e-11,15,5,2,8,3.6,10\n"
                     for i,p in enumerate(self.paths, 1)]

    def test_paired_valid(self):
        got = compare_naca_amg_v5.read_certified_replay(
            self.header + "".join(self.rows), self.paths, 1e-10)
        self.assertEqual(got["iterations"], 27)

    def test_paired_rejects_wrong_capture_and_order(self):
        broken = self.rows.copy()
        broken[0] = broken[0].replace("capture_1", "capture_8")
        with self.assertRaises(ValueError):
            compare_naca_amg_v5.read_certified_replay(
                self.header + "".join(broken), self.paths, 1e-10)
        broken = self.rows[::-1]
        with self.assertRaises(ValueError):
            compare_naca_amg_v5.read_certified_replay(
                self.header + "".join(broken), self.paths, 1e-10)

    def test_paired_rejects_false_residual(self):
        with self.assertRaises(ValueError):
            compare_naca_amg_v5.read_certified_replay(
                (self.header + "".join(self.rows)).replace("2e-11", "1e-7", 1),
                self.paths, 1e-10)

    def history(self):
        header = ("iteration,continuity_normalized,momentum_rhs_relative,"
                  "mass_normalized_imbalance,mass_local_linf,"
                  "flux_velocity_mismatch_normalized,pressure_linear_iterations,"
                  "momentum_linear_iterations,nonfinite_mass_faces,nonfinite_velocities\n")
        return header + "1,1e-5,2e-3,3e-4,4e-4,5e-4,9,3,0,0\n"

    def test_physics_valid(self):
        baseline = self.history()
        self.assertEqual(max(compare_naca_outer_histories.compare_histories(
            baseline,baseline,expected_iterations=1).values()), 0.0)

    def test_physics_rejects_mismatch_and_nonfinite(self):
        baseline = self.history()
        for altered in (baseline.replace("2e-3", "2e-2"),
                        baseline.replace("1e-5", "nan"),
                        baseline.replace(",9,3,", ",8,3,")):
            with self.subTest(altered=altered):
                with self.assertRaises(ValueError):
                    compare_naca_outer_histories.compare_histories(
                        baseline, altered, expected_iterations=1)


if __name__ == "__main__":
    unittest.main()
