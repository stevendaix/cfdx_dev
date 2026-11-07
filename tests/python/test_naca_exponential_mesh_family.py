"""Opt-in, physically two-dimensional NACA O-grid family regression."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'scripts' / 'generate_2d_validation_meshes.py'
SPEC = importlib.util.spec_from_file_location('naca_exponential_generator', SOURCE)
assert SPEC and SPEC.loader
MESHGEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MESHGEN)
sys.path.insert(0, str(ROOT / 'scripts'))
from naca_mesh_similarity_gate import verify_mesh_family  # noqa: E402


class ExponentialGridFamilyTests(unittest.TestCase):
    def test_exponential_family_refines_tangential_and_first_layer_together(self):
        rows = [MESHGEN.naca_o_grid_quality(n, nr, radial_stretch=6.3)
                for n, nr in ((32, 12), (64, 24), (128, 48))]
        report = verify_mesh_family(rows)
        self.assertEqual(report['status'], 'COMPATIBLE_WITH_SINGLE_RATIO_GCI', report)
        self.assertTrue(all(r['inverted_or_degenerate_cells'] == 0 for r in rows))
        self.assertLess(max(r['maximum_nonorthogonality_degrees'] for r in rows), 65.0)
        self.assertTrue(all(1.0 < p['first_layer_ratio'] < 2.5 for p in report['pairs']))

    def test_stretching_differences_block_a_false_gci(self):
        rows = [dict(n_surface=n, n_radial=nr, radius=20.0, radial_power=1.8,
                     radial_stretch=6.3,
                     minimum_first_layer_spacing=delta, inverted_or_degenerate_cells=0)
                for n, nr, delta in ((64,24,0.01),(128,48,0.005),(256,96,0.0025))]
        self.assertEqual(verify_mesh_family(rows)['status'], 'COMPATIBLE_WITH_SINGLE_RATIO_GCI')
        rows[1]['radial_stretch'] = 5.5
        self.assertIn('different_geometry_or_grading', verify_mesh_family(rows)['reasons'])

    def test_old_default_and_ambiguous_parameter_guard(self):
        self.assertIn('radial_stretch', MESHGEN.naca_o_grid_quality(16, 8))
        with self.assertRaisesRegex(ValueError, 'cannot be combined'):
            MESHGEN.naca_o_grid(16, 8, radial_power=2.0, radial_stretch=6.3)
        with self.assertRaisesRegex(ValueError, 'radial_stretch'):
            MESHGEN.naca_o_grid(16, 8, radial_stretch=float('nan'))


if __name__ == '__main__':
    unittest.main()
