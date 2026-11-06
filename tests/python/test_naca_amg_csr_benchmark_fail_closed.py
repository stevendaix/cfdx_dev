"""Fail-closed checks for the paired NACA V3/V4 benchmark parser."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / "tools/benchmarks/compare_naca_amg_csr.py"
spec = importlib.util.spec_from_file_location("compare_naca_amg_csr", SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestReplayCertification(unittest.TestCase):
    def write_runner(self, location: Path, status: str, residual: str):
        output = location / "fake_replay.py"
        output.write_text(
            "#!/usr/bin/env python3\n"
            "print('status,true_relative_residual')\n"
            f"print('{status},{residual}')\n", encoding="utf-8"
        )
        output.chmod(0o755)
        return output

    def test_explicit_failure_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            exe = self.write_runner(Path(d), "FAIL", "1e-12")
            with self.assertRaisesRegex(RuntimeError, "uncertified"):
                module.run(exe, "fgmres_boomer_fast", [Path("fake.manifest")], 5.0)

    def test_false_convergence_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            exe = self.write_runner(Path(d), "PASS", "2e-9")
            with self.assertRaisesRegex(RuntimeError, "uncertified"):
                module.run(exe, "fgmres_boomer_fast", [Path("fake.manifest")], 5.0)

    def test_certified_row_is_accepted(self):
        with tempfile.TemporaryDirectory() as d:
            exe = self.write_runner(Path(d), "PASS", "9e-12")
            rows = module.run(exe, "fgmres_boomer_fast", [Path("fake.manifest")], 5.0)
            self.assertEqual(len(rows), 1)


if __name__ == "__main__":
    unittest.main()
