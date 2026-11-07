"""Fail-closed contract tests for experimental coarse-to-fine initial guesses."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/naca_prolong_initial_state.py"
SPEC = importlib.util.spec_from_file_location("naca_prolong", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class ProlongTests(unittest.TestCase):
    def setUp(self):
        coords = np.array([[x, y, 0.0] for x in (0, 0.5, 1) for y in (0, 0.5, 1)])
        self.coarse = np.column_stack((coords, np.full(len(coords), 1.3),
                                       np.full(len(coords), -0.4), np.full(len(coords), 0.2)))
        self.fine = np.array([[0.25, 0.25, 0], [0.8, 0.4, 0]])

    def test_constant_solution_preserved(self):
        result = MODULE.prolong(self.coarse, self.fine)
        np.testing.assert_allclose(result[:, 3:], np.tile([1.3, -0.4, 0.2], (len(result), 1)), atol=1e-15)
        np.testing.assert_allclose(result[:, :3], self.fine)

    def test_exact_coincident_centers(self):
        self.coarse[0, 3:] = [2, 3, 4]
        result = MODULE.prolong(self.coarse, self.coarse[:1, :3])
        np.testing.assert_allclose(result[:, 3:], [[2, 3, 4]], atol=1e-14)

    def test_outside_domain_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.prolong(self.coarse, np.array([[99, 0, 0]]))

    def test_nan_rejected(self):
        self.coarse[0, 3] = np.nan
        with self.assertRaises(ValueError):
            MODULE.prolong(self.coarse, self.fine)

    def test_invalid_neighbours_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.prolong(self.coarse, self.fine, neighbors=2)

    def test_nested_affine_exact_on_curved_naca_grid(self):
        # The synthetic points follow an annular ring-major NACA-like grid.
        # Even a severe radial stretching must preserve every affine field.
        mc, mf = 9, 17
        nr_c, nr_f = 4, 8
        def points(nr: int, m: int):
            rows = []
            for j in range(nr):
                r = 1.0 + 0.25 * ((j + 0.5) / nr) ** 2
                for i in range(m):
                    theta = 2 * np.pi * (i + 0.5) / m
                    rows.append((r*np.cos(theta), r*np.sin(theta), 0.5))
            return np.asarray(rows)
        co = points(nr_c, mc)
        fi = points(nr_f, mf)
        coarse = np.column_stack((co, 0.2 + co[:, 0] - 2*co[:, 1],
                                  0.1 - 3*co[:, 0] + 0.5*co[:, 1],
                                  1.7 + 4*co[:, 1]))
        expected = np.column_stack((0.2 + fi[:, 0] - 2*fi[:, 1],
                                    0.1 - 3*fi[:, 0] + 0.5*fi[:, 1],
                                    1.7 + 4*fi[:, 1]))
        result = MODULE.prolong_nested_o_grid(
            coarse, fi, coarse_surface_n=4, fine_surface_n=8)
        np.testing.assert_allclose(result[:, :3], fi, atol=0, rtol=0)
        np.testing.assert_allclose(result[:, 3:], expected, atol=3e-12, rtol=0)

    def test_nested_rejects_false_refinement_or_broken_stencil(self):
        co = np.zeros((9*4, 6))
        fi = np.zeros((17*8, 3))
        with self.assertRaisesRegex(ValueError, "surface count must double"):
            MODULE.prolong_nested_o_grid(
                co, fi, coarse_surface_n=4, fine_surface_n=7)
        with self.assertRaisesRegex(ValueError, "degenerate coarse stencil"):
            MODULE.prolong_nested_o_grid(
                co, fi, coarse_surface_n=4, fine_surface_n=8)

    def test_nested_rejects_unrelated_or_nonplanar_grids(self):
        # Radial sizes match, but unrelated geometry must not be copied.
        mc, mf = 9, 17
        cr = np.array([[np.cos(2*np.pi*i/mc), np.sin(2*np.pi*i/mc), 0.0]
                       for _ in range(4) for i in range(mc)])
        fi = np.array([[np.cos(2*np.pi*i/mf)+100,
                        np.sin(2*np.pi*i/mf), 0.0]
                       for _ in range(8) for i in range(mf)])
        with self.assertRaisesRegex(ValueError, "unrelated mesh geometry"):
            MODULE.prolong_nested_o_grid(
                np.column_stack((cr, np.ones((len(cr), 3)))), fi,
                coarse_surface_n=4, fine_surface_n=8)
        fi[:, 2] = np.arange(len(fi))
        with self.assertRaisesRegex(ValueError, "planar extrusions"):
            MODULE.prolong_nested_o_grid(
                np.column_stack((cr, np.ones((len(cr), 3)))), fi,
                coarse_surface_n=4, fine_surface_n=8)

    def test_input_file_never_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coarse_path = root / "coarse.csv"
            fine_path = root / "fine.csv"
            with coarse_path.open("w") as stream:
                stream.write("cell,x,y,z,ux,uy,p\n")
                for i, row in enumerate(self.coarse):
                    stream.write(f"{i},"+",".join(map(str,row))+"\n")
            with fine_path.open("w") as stream:
                stream.write("cell,x,y,z\n")
                for i, row in enumerate(self.fine):
                    stream.write(f"{i},"+",".join(map(str,row))+"\n")
            original = coarse_path.read_bytes()
            with self.assertRaisesRegex(ValueError, "differ from input"):
                MODULE.create(coarse_path, fine_path, coarse_path)
            self.assertEqual(coarse_path.read_bytes(), original)

    def test_bad_csv_schema_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "data.csv"
            p.write_text("cell,x,y\n0,0,0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                MODULE.read_csv(p, ("cell", "x", "y", "z"))

    def test_missing_cell_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "data.csv"
            p.write_text("cell,x,y,z\n1,0,0,0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                MODULE.read_csv(p, ("cell", "x", "y", "z"))


if __name__ == "__main__":
    unittest.main()
