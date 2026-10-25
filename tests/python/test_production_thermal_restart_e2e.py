"""Thermal state must survive the real production executable's checkpoint."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

from test_solid_conduction_contract_e2e import write_cube


class ProductionThermalRestart(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.solver = Path(os.environ["CFDX_PRODUCTION_SOLVER"])
        if not cls.solver.is_file():
            raise RuntimeError(f"production executable is missing: {cls.solver}")

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cfdx_temperature_restart_")
        self.root = Path(self.tmp.name)
        self.case = self.root / "case.cfdx.h5"

        def configure(setup):
            setup.update(energy_model="conjugate", turbulence_model="laminar")
            setup["materials"][0].update(
                density=1.0, dynamic_viscosity=0.1, specific_heat=1.0
            )
            setup["numerics"] = {
                "selection": {
                    "entries": [
                        {
                            "family": "gradient",
                            "configuration_key": "numerics.gradient.gauss",
                        },
                        {
                            "family": "convection",
                            "configuration_key": "numerics.convection.upwind",
                        },
                        {
                            "family": "pressure_velocity",
                            "configuration_key": "pressure_velocity.simple",
                        },
                        {
                            "family": "linear_solver",
                            "configuration_key": "linear.fgmres",
                        },
                    ],
                    "required_families": [
                        "gradient",
                        "convection",
                        "pressure_velocity",
                    ],
                }
            }

        write_cube(self.case, configure)

    def tearDown(self):
        self.tmp.cleanup()

    def run_solver(self, name, restart=None):
        output = self.root / name
        command = [
            str(self.solver),
            "--mesh",
            str(self.case),
            "--output-dir",
            str(output),
            "--iterations",
            "4",
        ]
        if restart is not None:
            command += ["--restart", str(restart)]
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=30, check=False
        )
        return result, output

    def test_checkpoint_preserves_temperature_and_warm_energy_initialization(self):
        first, first_dir = self.run_solver("cold")
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        checkpoint = first_dir / "restart.dat"
        text = checkpoint.read_text()
        self.assertIn("field T 1\n", text)
        self.assertEqual(float(text.split("field T 1\n", 1)[1].split()[0]), 350.0)
        # The checkpoint is the thermal initializer. A cold solve would reject
        # this absent positive reference; a restart must consume its saved T.
        import h5py

        with h5py.File(self.case, "r+") as h5:
            setup = json.loads(h5.attrs["case_setup_json"])
            setup["ref_temperature"] = 0.0
            h5.attrs["case_setup_json"] = json.dumps(setup)
        second, second_dir = self.run_solver("warm", checkpoint)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        warm_energy = json.loads((second_dir / "energy_convergence.json").read_text())
        self.assertTrue(warm_energy["converged"])
        self.assertTrue(warm_energy["temperature_admissibility_passed"])
        self.assertEqual(warm_energy["nonpositive_temperature"], 0)
        self.assertEqual(warm_energy["iterations"], 1)
        restart_text = (second_dir / "restart.dat").read_text()
        self.assertEqual(
            float(restart_text.split("field T 1\n", 1)[1].split()[0]), 350.0
        )

    def test_thermal_restart_requires_temperature_instead_of_cold_fallback(self):
        first, output = self.run_solver("cold")
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        checkpoint = self.root / "missing_temperature.dat"
        text = (output / "restart.dat").read_text()
        checkpoint.write_text(text.split("field T 1\n", 1)[0])
        result, rejected_dir = self.run_solver("rejected", checkpoint)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("requested T but DAT has none", result.stderr)
        self.assertFalse((rejected_dir / "restart.dat").exists())
        execution = json.loads((rejected_dir / "execution.json").read_text())
        self.assertFalse(execution["converged"])

    def test_negative_restart_temperature_is_rejected_before_flow(self):
        first, output = self.run_solver("cold")
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        checkpoint = self.root / "negative_temperature.dat"
        text = (output / "restart.dat").read_text()
        checkpoint.write_text(text.split("field T 1\n", 1)[0] + "field T 1\n-1\n")
        result, rejected_dir = self.run_solver("rejected", checkpoint)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("restarted temperature", result.stderr)
        self.assertFalse((rejected_dir / "restart.dat").exists())

    def test_physical_transient_flag_cannot_be_solved_as_steady(self):
        import h5py

        with h5py.File(self.case, "r+") as h5:
            setup = json.loads(h5.attrs["case_setup_json"])
            setup.update(transient=True, time_step=0.1, end_time=1.0)
            h5.attrs["case_setup_json"] = json.dumps(setup)
        result, output = self.run_solver("rejected_time")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("transient evolution", result.stderr)
        self.assertFalse((output / "restart.dat").exists())

    def test_energy_exception_cannot_report_successful_carrier_as_overall_convergence(
        self,
    ):
        import h5py

        for missing in ("ref_temperature", "specific_heat"):
            with self.subTest(missing=missing):
                with h5py.File(self.case, "r+") as h5:
                    setup = json.loads(h5.attrs["case_setup_json"])
                    setup["ref_temperature"] = 300.0
                    setup["materials"][0]["specific_heat"] = 1.0
                    if missing == "ref_temperature":
                        setup[missing] = 0.0
                    else:
                        setup["materials"][0].pop(missing)
                    h5.attrs["case_setup_json"] = json.dumps(setup)
                result, output = self.run_solver("missing_" + missing)
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                execution = json.loads((output / "execution.json").read_text())
                self.assertFalse(execution["converged"])
                self.assertEqual(execution["process_exit_code"], 2)
                self.assertIn(execution["convergence_reason"], result.stderr)

    def test_conservative_negative_temperature_is_not_physically_accepted(self):
        import h5py

        with h5py.File(self.case, "r+") as h5:
            setup = json.loads(h5.attrs["case_setup_json"])
            # Two unit-area walls at 300/400 K, lambda=2 and d=0.5:
            # 8*T = 2800 + Su. Su=-3600 yields exactly T=-100 K.
            setup["materials"][0]["extra_properties"] = {"volumetric_power": "-3600"}
            h5.attrs["case_setup_json"] = json.dumps(setup)
        result, output = self.run_solver("negative_final")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        energy = json.loads((output / "energy_convergence.json").read_text())
        summary = json.loads((output / "energy_summary.json").read_text())
        execution = json.loads((output / "execution.json").read_text())
        self.assertTrue(energy["linear_solver_converged"])
        self.assertTrue(energy["global_energy_balance_passed"])
        self.assertTrue(energy["local_energy_balance_passed"])
        self.assertFalse(energy["temperature_admissibility_passed"])
        self.assertEqual(energy["nonpositive_temperature"], 1)
        self.assertEqual(summary["min_temperature"], -100.0)
        self.assertFalse(energy["converged"])
        self.assertFalse(execution["converged"])
        self.assertIn("temperature admissibility", execution["convergence_reason"])


if __name__ == "__main__":
    unittest.main(verbosity=2)

