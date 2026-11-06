"""An incomplete NACA solve must not leave old forces labelled converged."""
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


def ctest_arguments():
    """Return the (binary, mesh) paths provided by the CTest registration.

    CTest runs this module as a script with the NACA qualification binary
    and the quick n64 mesh as command-line arguments. The full-tree pytest
    regression collects the module without those arguments, so pytest's
    own options must never be misread as the solver path: without the
    CTest-provided binary and mesh the fail-closed gate skips instead of
    executing a non-existent file.
    """
    if len(sys.argv) >= 3:
        binary = Path(sys.argv[1])
        mesh = Path(sys.argv[2])
        if binary.is_file() and mesh.is_file():
            return binary, mesh
    raise unittest.SkipTest(
        "certificate gate requires the CTest-provided binary and mesh"
    )


class CertificateGate(unittest.TestCase):
    def test_failure_records_history_and_removes_stale_force_result(self):
        binary, mesh = ctest_arguments()
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            forces = directory / "forces.csv"
            history = directory / "physical_history.csv"
            forces.write_text("stale converged forces must be removed")
            env = os.environ.copy()
            env["CFDX_NACA_PHYSICAL_HISTORY_CSV"] = str(history)
            env["CFDX_NACA_CONVERGED_FORCE_CSV"] = str(forces)
            env["CFDX_NACA_PRESSURE_SOLVER"] = "fgmres_boomer"
            proc = subprocess.run(
                [str(binary), "--single", str(mesh), "0.6", "0.25", "1"],
                env=env, capture_output=True, text=True, timeout=40,
                check=False,
            )
            self.assertEqual(proc.returncode, 1, proc.stderr[-1500:])
            self.assertIn("NACA0012_QUALIFICATION: FAIL", proc.stderr)
            self.assertFalse(forces.exists(), "stale force certificate survived a failure")
            self.assertTrue(history.is_file())
            with history.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 1)
            self.assertEqual(int(rows[0]["iteration"]), 1)
            self.assertNotIn("NACA_CONVERGED_FORCE_CERTIFICATE", proc.stdout)

    def test_invalid_mesh_clears_stale_force_and_state_certificates(self):
        binary, _ = ctest_arguments()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            force = path / "old_forces.csv"
            state = path / "old_solution.csv"
            profile = path / "old_face_profile.csv"
            history = path / "old_physical_history.csv"
            force.write_text("old_force_success")
            state.write_text("old_state_success")
            profile.write_text("old_profile_success")
            history.write_text("old_history_success")
            env = os.environ.copy()
            env["CFDX_NACA_CONVERGED_FORCE_CSV"] = str(force)
            env["CFDX_NACA_EXPORT_CONVERGED_STATE_CSV"] = str(state)
            env["CFDX_NACA_FACE_PROFILES_CSV"] = str(profile)
            env["CFDX_NACA_PHYSICAL_HISTORY_CSV"] = str(history)
            proc = subprocess.run(
                [str(binary), "--single", str(path / "nonexistent.h5"),
                 "0.6", "0.25", "1"],
                env=env, capture_output=True, text=True, timeout=15,
                check=False,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertFalse(force.exists(), "stale force survived mesh loading failure")
            self.assertFalse(state.exists(), "stale state survived mesh loading failure")
            self.assertFalse(profile.exists(), "stale wall profile survived mesh loading failure")
            self.assertFalse(history.exists(), "stale history survived mesh loading failure")

    def test_output_preflight_protects_restart_input(self):
        binary, _ = ctest_arguments()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            restart = path / "restart.csv"
            restart.write_text("critical original input")
            env = os.environ.copy()
            env["CFDX_NACA_IMPORT_INITIAL_STATE_CSV"] = str(restart)
            env["CFDX_NACA_CONVERGED_FORCE_CSV"] = str(restart)
            proc = subprocess.run(
                [str(binary), "--single", str(path / "missing.h5"),
                 "0.6", "0.25", "1"],
                env=env, capture_output=True, text=True, timeout=15,
                check=False,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("conflicting NACA output path", proc.stderr)
            self.assertEqual(restart.read_text(), "critical original input")


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]], verbosity=2)
