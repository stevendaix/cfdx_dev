"""Adversarial checks: no cap is qualified on malformed or false-positive CSV."""
import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / "tools/benchmarks/qualify_naca_amg_cap_sweep.py"
spec = importlib.util.spec_from_file_location("amg_cap_tuner", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestFailClosed(unittest.TestCase):
    def setUp(self):
        self.paths = [Path(f"/tmp/amg_capture_{i}.manifest") for i in range(1, 4)]
        head = "manifest,status,sequence,rows,iterations,true_relative_residual,total_ms,full_setup_ms,numeric_update_ms,krylov_ms,amg_operator_complexity,amg_levels\n"
        self.lines = [
            f"{name},PASS,{i},12336,12,8e-11,20.0,7.0,3.0,10.0,3.6,10\n"
            for i, name in enumerate(self.paths, start=1)
        ]
        self.csv = head + "".join(self.lines)

    def test_valid_three_corrections(self):
        metrics = module.read_replay(self.csv, self.paths, 1e-10)
        self.assertEqual(metrics["iterations"], 36)
        self.assertEqual(metrics["total_ms"], 60.0)

    def test_not_all_corrections_reported(self):
        with self.assertRaises(ValueError):
            module.read_replay(self.csv.rsplit("\n", 2)[0] + "\n", self.paths, 1e-10)

    def test_false_pass_and_infinite_residual_rejected(self):
        for altered in (self.csv.replace("8e-11", "4e-9", 1),
                        self.csv.replace("8e-11", "nan", 1),
                        self.csv.replace("PASS,2", "FAIL,2")):
            with self.subTest(altered=altered[:90]):
                with self.assertRaises(ValueError):
                    module.read_replay(altered, self.paths, 1e-10)

    def test_capture_identity_and_order_enforced(self):
        bad = self.csv.replace("amg_capture_2.manifest", "amg_capture_9.manifest")
        with self.assertRaises(ValueError):
            module.read_replay(bad, self.paths, 1e-10)

    def test_nonphysical_benchmarks_rejected(self):
        for altered in (self.csv.replace("20.0,7.0", "-1,7.0"),
                        self.csv.replace("3.6,10", "nan,10")):
            with self.assertRaises(ValueError):
                module.read_replay(altered, self.paths, 1e-10)


if __name__ == "__main__":
    unittest.main()
