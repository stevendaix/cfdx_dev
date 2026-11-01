"""No fake three-grid NACA wall-shear qualification."""

import importlib.util
import math
from pathlib import Path
import unittest


PATH = Path(__file__).resolve().parents[2] / "scripts/naca_wall_shear_independent.py"
SPEC = importlib.util.spec_from_file_location("wall_shear_independent", PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def good_input():
    return {
        "case": "NACA0012_RE1000_A0",
        "results": [
            {
                "N": n,
                "status": "PASS",
                "converged": True,
                "continuity": 1e-9,
                "wall_faces": 2*n + 1,
                "covered_wall_faces": 2*n + 1,
                "cd_viscous_face": v,
                "cd_viscous_quadratic": w,
            }
            for n, v, w in (
                (64, 0.085, 0.086),
                (128, 0.080, 0.0806),
                (256, 0.07875, 0.07925),
            )
        ],
    }


class WallShearGateTests(unittest.TestCase):
    def test_three_grid_coherent(self):
        result = MODULE.verify(good_input())
        self.assertEqual(result["status"], "GRID_AND_METHOD_CONSISTENT")
        self.assertEqual(result["reference_validation"], "NOT_PERFORMED")
        self.assertTrue(math.isfinite(result["studies"]["cd_viscous_face"]["observed_order"]))

    def test_missing_n128(self):
        sample = good_input()
        sample["results"].pop(1)
        with self.assertRaises(ValueError):
            MODULE.verify(sample)

    def test_unconverged_n256(self):
        sample = good_input()
        sample["results"][2]["status"] = "TIMEOUT"
        with self.assertRaises(ValueError):
            MODULE.verify(sample)

    def test_false_convergence(self):
        sample = good_input()
        sample["results"][2]["converged"] = False
        with self.assertRaises(ValueError):
            MODULE.verify(sample)

    def test_missing_wall_faces(self):
        sample = good_input()
        sample["results"][1]["covered_wall_faces"] -= 1
        with self.assertRaises(ValueError):
            MODULE.verify(sample)

    def test_nonfinite_force(self):
        sample = good_input()
        sample["results"][2]["cd_viscous_quadratic"] = float("nan")
        with self.assertRaises(ValueError):
            MODULE.verify(sample)

    def test_methods_disagree(self):
        sample = good_input()
        sample["results"][2]["cd_viscous_quadratic"] = 0.12
        with self.assertRaises(ValueError):
            MODULE.verify(sample)

    def test_nonmonotone_refinement(self):
        sample = good_input()
        sample["results"][2]["cd_viscous_face"] = 0.081
        with self.assertRaises(ValueError):
            MODULE.verify(sample)

    def test_missing_reference_never_becomes_external_validation(self):
        result = MODULE.verify(good_input())
        self.assertNotEqual(result["reference_validation"], "PASS")


if __name__ == "__main__":
    unittest.main()
