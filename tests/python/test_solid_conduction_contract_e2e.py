"""Exercise physics selection through the real solid-conduction executable."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import h5py
import numpy as np


def write_cube(path: Path, change=None) -> None:
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
    faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
             [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]]
    names = ["bottom", "top", "front", "back", "left", "right"]
    setup = {
        "energy_model": "solid_conduction", "solver_mode": "steady",
        "ref_temperature": 300.0,
        "materials": [{"name": "solid", "thermal_conductivity": 2.0}],
        "boundary_conditions": [
            {"patch_name": name, "type": "wall", "value_type": "wall_thermal",
             "wall_temperature": 300.0 if name == "left" else 400.0 if name == "right" else 0.0}
            for name in names],
    }
    if change is not None:
        change(setup)
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["format_version"] = "1"
        h5.attrs["schema_version"] = "1"
        h5.attrs["mesh_topology"] = "cfdx-csr-v1"
        h5.attrs["case_setup_json"] = json.dumps(setup)
        h5.attrs["boundary_patches"] = ";".join(
            f"{name}:{index}:1:0" for index, name in enumerate(names))
        h5.create_dataset("points", data=points)
        h5.create_dataset("face_vertices", data=np.array(faces, dtype=np.uint64).ravel())
        h5.create_dataset("face_offsets", data=np.arange(0, 25, 4, dtype=np.uint64))
        h5.create_dataset("owner", data=np.zeros(6, dtype=np.uint64))
        h5.create_dataset("neighbour", data=np.full(6, -1, dtype=np.int64))
        h5.create_dataset("cell_faces", data=np.arange(6, dtype=np.uint64))
        h5.create_dataset("cell_offsets", data=np.array([0, 6], dtype=np.uint64))
        h5.create_dataset("patch_face_ids", data=np.arange(6, dtype=np.uint64))
        h5.create_dataset("patch_face_offsets", data=np.arange(7, dtype=np.uint64))


class SolidConductionApplicationContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.solver = Path(os.environ["CFDX_THERMAL_SOLVER"])
        if not cls.solver.is_file():
            raise RuntimeError(f"thermal executable is missing: {cls.solver}")

    def run_case(self, change=None):
        with tempfile.TemporaryDirectory(prefix="cfdx_solid_contract_") as work:
            root = Path(work)
            write_cube(root / "case.cfdx.h5", change)
            result = subprocess.run(
                [str(self.solver), "--mesh", str(root / "case.cfdx.h5"),
                 "--output-dir", str(root / "output")],
                capture_output=True, text=True, timeout=30, check=False)
            execution = json.loads((root / "output/execution.json").read_text())
            summary_path = root / "output/thermal_summary.json"
            summary = json.loads(summary_path.read_text()) if summary_path.exists() else None
            return result, execution, summary

    def assert_rejected(self, change, message):
        result, execution, summary = self.run_case(change)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(execution["process_exit_code"], 2)
        self.assertFalse(execution["converged"])
        self.assertIsNone(summary)
        self.assertIn(message, result.stderr)

    def test_supported_solid_case_solves(self):
        result, execution, summary = self.run_case()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(execution["converged"])
        self.assertIsNotNone(summary)
        self.assertIn("Temperature min 350 max 350", result.stdout)

    def test_radiation_is_not_silently_ignored(self):
        self.assert_rejected(lambda case: case.update(radiation_model="S2S"), "radiation")

    def test_transient_evolution_is_not_silently_ignored(self):
        self.assert_rejected(lambda case: case.update(transient=True, time_step=0.1), "transient")

    def test_turbulence_is_not_silently_ignored(self):
        self.assert_rejected(lambda case: case.update(turbulence_model="k_omega_sst"), "turbulence")

    def test_canonical_temporal_selection_is_not_silently_ignored(self):
        self.assert_rejected(lambda case: case.update(numerics={"selection": {
            "entries": [{"family": "temporal", "configuration_key": "numerics.temporal.bdf2"}],
            "required_families": ["temporal"],
        }}), "transient evolution")

    def test_multiphase_is_not_silently_ignored(self):
        self.assert_rejected(lambda case: case.update(multiphase_model="VOF"), "multiphase")

    def test_fluid_inlet_is_not_replaced_by_adiabatic_wall(self):
        def change(case):
            case["boundary_conditions"][0].update(type="inlet", velocity_magnitude=1.0)
        self.assert_rejected(change, "boundary physics")

    def test_absent_conductivity_is_not_replaced_by_air(self):
        self.assert_rejected(lambda case: case["materials"][0].pop("thermal_conductivity"),
                             "conductivity")


if __name__ == "__main__":
    unittest.main(verbosity=2)
