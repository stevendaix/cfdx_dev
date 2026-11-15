"""Tests for native CFDX cell-field maps and honest VTK/PyVista export."""

from __future__ import annotations

import importlib.util
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "plot_bfs_v83_pyvista.py"
SPEC = importlib.util.spec_from_file_location("plot_bfs_v83_pyvista", SCRIPT)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class TestLimits(unittest.TestCase):
    def fixtures(self):
        base = {k: np.asarray([0.0, 0.5], dtype=float) for k in M.FIELDS}
        base["u_x"][:] = [36.0, 50.0]
        base["u_mag"][:] = [40.0, 51.0]
        base["p_kpa"][:] = [-25.0, -0.1]
        base["omega_log"][:] = [2.3, 7.0]
        base["k"][:] = [1e-4, 2e-3]
        return base

    def test_limits_center_reference_zero(self):
        sc = M.color_limits([self.fixtures()])
        self.assertEqual(sc["u_y"], (-13.0, 13.0))
        self.assertEqual(sc["p_kpa"], (-30.0, 0.0))

    def test_backflow_does_not_get_clipped(self):
        base = self.fixtures()
        base["u_x"][0] = -18.0
        sc = M.color_limits([base])
        self.assertLessEqual(sc["u_x"][0], -18.0)

    def test_multiple_states_share_scales(self):
        a = self.fixtures()
        b = self.fixtures()
        b["u_mag"][0] = 60.0
        sc = M.color_limits([a, b])
        self.assertEqual(sc["u_mag"][1], 60.0)

    def test_nonfinite_rejected_from_vtu_and_invalid_height(self):
        with self.assertRaises(ValueError):
            M.load_vtu(Path("unused"), h=0)

    def test_frame_geometry_aspect(self):
        bounds = M._frame((-5, 25, 0, 4), 1800, 360)
        self.assertAlmostEqual((bounds[1] - bounds[0]) / (bounds[3] - bounds[2]), 5)
        self.assertLessEqual(bounds[0], -5)
        self.assertGreaterEqual(bounds[1], 25)


class TestNativeCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = os.environ.get("CFDX_V83_NATIVE_VTU", "")
        cls.native = Path(root) if root else None

    def test_native_cell_fields_and_positive_speed(self):
        if self.native is None:
            self.skipTest("Set CFDX_V83_NATIVE_VTU for real mesh test")
        mesh, arrays = M.load_vtu(self.native)
        self.assertEqual(mesh.GetNumberOfCells(), 20540)
        self.assertTrue(np.array_equal(arrays["p_kpa"], arrays["p"] / 1000))
        self.assertTrue(np.all(arrays["u_mag"] >= 0))
        self.assertEqual(int(np.sum(arrays["u_x"] < 0)), 0)
        self.assertEqual(mesh.GetCellData().GetArray("u_x").GetNumberOfTuples(), 20540)

    def test_native_pyvista_source_is_present(self):
        self.assertIn("import pyvista as pv", SCRIPT.read_text())
        self.assertIn('preference="cell"', SCRIPT.read_text())

    def test_native_render_writes_real_map(self):
        if self.native is None:
            self.skipTest("Set CFDX_V83_NATIVE_VTU for real map test")
        mesh, a = M.load_vtu(self.native)
        surface = M._vtk_surface(mesh)
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "render.png"
            extent = M.render_vtk(
                surface,
                "p_kpa",
                (-30, 0),
                "coolwarm",
                (-2, 10, 0, 3),
                p,
                window=(720, 240),
            )
            self.assertTrue(p.is_file())
            self.assertGreater(p.stat().st_size, 1000)
            self.assertLessEqual(extent[0], -2.0)


if __name__ == "__main__":
    unittest.main()
