"""Fail-closed contract tests for scripts/naca_controlled_campaign.py."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "naca_controlled_campaign.py"
SPEC = importlib.util.spec_from_file_location("naca_controlled_campaign", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
CAMPAIGN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CAMPAIGN)

parse_result = CAMPAIGN.parse_result
run_campaign = CAMPAIGN.run_campaign


class ControlledNacaTests(unittest.TestCase):
    sample = (
        "NACA0012_RE1000_A0 level=0 Cd=0.09 Cl=0 "
        "Cm=0 Cd_pressure=0.03 Cd_viscous=0.06 "
        "continuity=1e-10 iterations=456\n"
        "NACA0012_SINGLE: PASS (one level only; not a mesh qualification)\n"
    )

    def test_valid_sample(self):
        self.assertAlmostEqual(parse_result(self.sample)["Cd"], 0.09)

    def test_missing_success_rejected(self):
        with self.assertRaises(ValueError):
            parse_result(self.sample.replace("NACA0012_SINGLE: PASS", ""))

    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            parse_result(self.sample.replace("Cd=0.09", "Cd=nan"))

    def test_timeout_is_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mesh = root / "mesh.h5"
            mesh.write_text("fixture")
            script = root / "slow.py"
            script.write_text(
                "#!/usr/bin/env python3\nimport time\ntime.sleep(20)\n"
            )
            script.chmod(0o755)
            result = run_campaign(script, [mesh, mesh, mesh], 0.1, root / "logs")
            self.assertEqual(result["three_grid_solver_gate"], "FAIL")
            self.assertEqual(result["results"][0]["reason"], "timeout")

    def test_nonzero_exit_is_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mesh = root / "mesh.h5"
            mesh.write_text("fixture")
            script = root / "fail.py"
            script.write_text(
                "#!/usr/bin/env python3\nimport sys\n"
                'print("NACA0012_SINGLE: PASS")\nsys.exit(4)\n'
            )
            script.chmod(0o755)
            result = run_campaign(script, [mesh, mesh, mesh], 2, root / "logs")
            self.assertEqual(result["three_grid_solver_gate"], "FAIL")
            self.assertEqual(result["results"][0]["exit_code"], 4)


if __name__ == "__main__":
    unittest.main()
