import importlib.util
from pathlib import Path
import tempfile
import unittest

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "naca_ab_qualification.py"
_spec = importlib.util.spec_from_file_location("cfdx_naca_ab", _SCRIPT)
assert _spec is not None and _spec.loader is not None
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

_PASS = (
    "NACA0012_RE1000_A0 level=0 Cd=0.091 Cl=0.0 Cm=0 "
    "Cd_pressure=0.06 Cd_viscous=0.031 continuity=1e-11 iterations=450\n"
    "NACA0012_SINGLE: PASS (one level only; not a mesh qualification)\n"
)


class TestNacaAB(unittest.TestCase):
    def test_true_success(self):
        self.assertEqual(_mod.parse_output(_PASS, 0)["iterations"], 450)

    def test_nonzero_exit_rejected(self):
        with self.assertRaises(ValueError):
            _mod.parse_output(_PASS, 1)

    def test_missing_marker_rejected(self):
        with self.assertRaises(ValueError):
            _mod.parse_output(_PASS.replace("NACA0012_SINGLE: PASS", "FAILED"), 0)

    def test_nan_rejected(self):
        with self.assertRaises(ValueError):
            _mod.parse_output(_PASS.replace("Cd=0.091", "Cd=nan"), 0)

    def test_nonconservative_force_rejected(self):
        with self.assertRaises(ValueError):
            _mod.parse_output(_PASS.replace("Cd_viscous=0.031", "Cd_viscous=0.001"), 0)

    def test_negative_continuity_rejected(self):
        with self.assertRaises(ValueError):
            _mod.parse_output(_PASS.replace("continuity=1e-11", "continuity=-1e-11"), 0)

    def test_zero_iteration_pass_marker_rejected(self):
        with self.assertRaises(ValueError):
            _mod.parse_output(_PASS.replace("iterations=450", "iterations=0"), 0)

    def test_no_false_pass_from_timeout(self):
        with tempfile.TemporaryDirectory() as td:
            dummy = Path(td) / "fake.py"
        with self.assertRaises(ValueError):
            _mod.parse_output(_PASS.replace("Cd_viscous=0.031", "Cd_viscous=0.001"), 0)

    def test_no_false_pass_from_timeout(self):
        with tempfile.TemporaryDirectory() as td:
            dummy = Path(td) / "fake.py"
            dummy.write_text("unused")
            # The fake executable emits a success marker but times out.
            # This must never count as a successful CFD qualification.
            exe = Path(td) / "fake.sh"
            exe.write_text("#!/bin/sh\necho 'NACA0012_SINGLE: PASS'\nsleep 2\n")
            exe.chmod(0o755)
            result = _mod.run_one(exe, dummy, "baseline", 1, Path(td) / "log.txt")
            self.assertEqual(result["status"], "FAIL")
            self.assertEqual(result["reason"], "TIMEOUT")


if __name__ == "__main__":
    unittest.main()
