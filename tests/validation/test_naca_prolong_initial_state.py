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
