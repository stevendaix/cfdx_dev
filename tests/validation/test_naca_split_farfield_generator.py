import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py


_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "generate_2d_validation_meshes.py"
_spec = importlib.util.spec_from_file_location("cfdx_naca_mesher", _SCRIPT)
assert _spec is not None and _spec.loader is not None
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


class TestFarfieldSplit(unittest.TestCase):
    def test_split_and_default_keep_all_faces(self):
        args = (16, 6, 20.0)
        _, _, baseline = _mod.naca_o_grid(*args)
        points, cells, candidate = _mod.naca_o_grid(*args, split_farfield=True)
        self.assertIn("farfield", baseline)
        self.assertNotIn("farfield", candidate)
        self.assertTrue(candidate["farfield_inlet"])
        self.assertTrue(candidate["farfield_outlet"])
        self.assertEqual(
            len(candidate["farfield_inlet"]) + len(candidate["farfield_outlet"]),
            len(baseline["farfield"]),
        )
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "candidate.h5"
            _mod.write_mesh(output, points, cells, candidate)
            with h5py.File(output, "r") as f:
                names = f.attrs["boundary_patches"]
                self.assertIn("farfield_inlet", names)
                self.assertIn("farfield_outlet", names)
                self.assertNotIn("farfield:", names)

    def test_invalid_geometry_fails_closed(self):
        with self.assertRaises(ValueError):
            _mod.naca_o_grid(16, 6, radius=0)


if __name__ == "__main__":
    unittest.main()
