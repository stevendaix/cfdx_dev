import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from generate_2d_validation_meshes import naca_o_grid  # noqa: E402
from naca_radial_sweep import sweep  # noqa: E402


class SweepTests(unittest.TestCase):
    def test_default_unchanged(self):
        self.assertEqual(
            naca_o_grid(16, 8),
            naca_o_grid(16, 8, radial_power=1.8),
        )

    def test_invalid_power_rejected(self):
        for power in (float("nan"), float("inf"), 0.99, 3.01):
            with self.subTest(power=power):
                with self.assertRaises(ValueError):
                    naca_o_grid(16, 8, radial_power=power)

    def test_sweep_and_metrics(self):
        result = sweep(16, 8, [1.5, 1.8])
        self.assertEqual(len(result["candidates"]), 2)
        self.assertAlmostEqual(
            result["candidates"][1]["first_layer_ratio_to_baseline"], 1.0
        )
        self.assertEqual(
            result["candidates"][1]["maximum_radial_growth_ratio"],
            result["baseline"]["maximum_radial_growth_ratio"],
        )
        self.assertEqual(result["candidates"][1]["radial_power"], 1.8)

    def test_invalid_sweep_threshold_rejected(self):
        with self.assertRaises(ValueError):
            sweep(16, 8, [1.8], max_nonorth=float("nan"))


if __name__ == "__main__":
    unittest.main()
