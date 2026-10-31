"""Fail-closed tests for the qualification matrix runner."""

import importlib.util
import sys
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts" / "run_qualification_matrix.py"
)
spec = importlib.util.spec_from_file_location("qualification_matrix", SCRIPT)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class TestMatrix(unittest.TestCase):
    def test_ctest_pass_requires_exact_test(self):
        log = (
            "1/1 Test #123: test_thermal_vv .... Passed 0.02 sec\n"
            "100% tests passed, 0 tests failed out of 1\n"
        )
        self.assertTrue(module.check_ctest_output(0, log, "test_thermal_vv"))
        self.assertFalse(module.check_ctest_output(0, log, "other"))
        self.assertFalse(module.check_ctest_output(1, log, "test_thermal_vv"))

    def test_ctest_fixtures_can_run_with_target(self):
        log = (
            "1/3 Test #7: mesh_fixture .... Passed 0.01 sec\n"
            "2/3 Test #8: test_thermal_vv .... Passed 0.02 sec\n"
            "3/3 Test #9: cleanup_fixture .... Passed 0.01 sec\n"
            "100% tests passed, 0 tests failed out of 3\n"
        )
        self.assertTrue(module.check_ctest_output(0, log, "test_thermal_vv"))
        self.assertFalse(module.check_ctest_output(0, log, "test_absent"))

    def test_empty_registry_cannot_pass(self):
        self.assertFalse(
            module.check_ctest_output(0, "No tests were found!!!", "anything")
        )
        self.assertFalse(
            module.check_ctest_output(
                0, "100% tests passed, 0 tests failed out of 0", "anything"
            )
        )

    def test_failed_and_notrun_cannot_pass(self):
        for line in ["***Failed", "***Not Run", "***Timeout"]:
            log = (
                f"1/1 Test #1: test_thermal_vv .... {line}\n"
                "0% tests passed, 1 tests failed out of 1\n"
            )
            self.assertFalse(
                module.check_ctest_output(0, log, "test_thermal_vv")
            )

    def test_long_gate_explicit(self):
        self.assertEqual(len(module.QUICK), 9)
        self.assertEqual(
            {c.name for c in module.LONG},
            {
                "test_naca0012_qualification",
                "test_naca0012_wall_shear_audit",
                "test_ghia_cavity_profile_study",
            },
        )
        self.assertTrue(all(c.level == "physical" for c in module.LONG))


if __name__ == "__main__":
    unittest.main()
