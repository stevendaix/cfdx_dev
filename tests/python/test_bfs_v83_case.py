import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError:
    h5py = None
from scripts.enable_bfs20540_v83_progress import REQUIRED, make_progress_case


@unittest.skipUnless(h5py is not None, "h5py required")
class BFS20540V83CaseTests(unittest.TestCase):
    def test_case_preserves_mesh_and_thresholds(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            s = root / "original.h5"
            out = root / "progress.h5"
            setup = {
                "transient": False,
                "numerics": {
                    "coupled_solver": "SIMPLE",
                    "raw_settings": dict(
                        REQUIRED,
                        pressure_linear_tolerance="1e-11",
                        sst_transport_max_iterations="10",
                    ),
                },
            }
            with h5py.File(s, "w") as h:
                h.attrs["n_cells"] = 20540
                h.attrs["source_case_name"] = "backwardFacingStep2D"
                h.attrs["mesh_hash"] = "same_mesh"
                h.attrs["topology_hash"] = "same_topology"
                h.attrs["case_setup_json"] = json.dumps(setup)
                h.create_dataset("cells", data=np.array([11, 22, 33], dtype="i8"))
            result = make_progress_case(s, out)
            self.assertFalse(result["physics_qualified"])
            with h5py.File(out, "r") as f:
                changed = json.loads(f.attrs["case_setup_json"])
                raw = changed["numerics"]["raw_settings"]
                self.assertEqual(raw["sst_outer_provisional"], "enabled")
                self.assertEqual(raw["v83_provisional_outer_merit_progress"], "enabled")
                self.assertEqual(raw["pressure_linear_tolerance"], "1e-11")
                self.assertEqual(f["cells"][:].tolist(), [11, 22, 33])
            with self.assertRaises(ValueError):
                make_progress_case(s, s)
            with h5py.File(s, "r+") as f:
                f.attrs["n_cells"] = 152
            with self.assertRaises(ValueError):
                make_progress_case(s, root / "wrong_mesh.h5")
            with h5py.File(s, "r+") as f:
                f.attrs["n_cells"] = 20540
                setup["numerics"]["raw_settings"]["v8_global_nonlinear_residual"] = (
                    "false"
                )
                f.attrs["case_setup_json"] = json.dumps(setup)
            with self.assertRaises(ValueError):
                make_progress_case(s, root / "missing_v8.h5")


if __name__ == "__main__":
    unittest.main()
