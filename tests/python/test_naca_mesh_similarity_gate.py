"""A GCI on one h must not hide disproportionate near-wall refinement."""
import importlib.util
import unittest
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[2] / "scripts/naca_mesh_similarity_gate.py"
SPEC = importlib.util.spec_from_file_location("naca_mesh_similarity_gate", SOURCE)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row(surface, radial, spacing):
    return {"n_surface": surface, "n_radial": radial, "radius": 20.0,
            "radial_power": 1.8, "minimum_first_layer_spacing": spacing,
            "inverted_or_degenerate_cells": 0}


class MeshFamilyTests(unittest.TestCase):
    def test_geometrically_similar_mesh_family_is_allowed(self):
        self.assertEqual(MODULE.verify_mesh_family([
            row(64, 24, 0.04), row(128, 48, 0.02),
            row(256, 96, 0.01)])["status"],
            "COMPATIBLE_WITH_SINGLE_RATIO_GCI")

    def test_wall_normal_refinement_mismatch_blocks_gci(self):
        report = MODULE.verify_mesh_family([
            row(64, 24, 0.05), row(128, 48, 0.0155),
            row(256, 96, 0.0047)])
        self.assertEqual(report["status"], "BLOCKED")
        self.assertIn("first_layer_refinement_not_similar", report["reasons"])

    def test_bad_geometry_and_inverted_cells_block(self):
        rows = [row(64,24,0.04),row(128,48,0.02),row(256,96,0.01)]
        rows[1]["radius"] = 10.0
        rows[2]["inverted_or_degenerate_cells"] = 1
        self.assertEqual(MODULE.verify_mesh_family(rows)["status"], "BLOCKED")

    def test_mixed_radial_distributions_are_not_comparable(self):
        rows = [row(64,24,0.04), row(128,48,0.02), row(256,96,0.01)]
        rows[0]["radial_refinement_strategy"] = "powerlaw_resampled"
        rows[1]["radial_refinement_strategy"] = "nested_base_24"
        rows[2]["radial_refinement_strategy"] = "nested_base_24"
        result = MODULE.verify_mesh_family(rows)
        self.assertEqual(result["status"], "BLOCKED")
        self.assertIn("inconsistent_radial_refinement_strategy", result["reasons"])

    def test_nonfinite_spacing_fails_closed(self):
        rows = [row(64,24,0.04),row(128,48,float("nan")),row(256,96,0.01)]
        with self.assertRaisesRegex(ValueError, "invalid mesh characteristic"):
        rows[2]["inverted_or_degenerate_cells"] = 1
        self.assertEqual(MODULE.verify_mesh_family(rows)["status"], "BLOCKED")

    def test_nonfinite_spacing_fails_closed(self):
        rows = [row(64,24,0.04),row(128,48,float("nan")),row(256,96,0.01)]
        with self.assertRaisesRegex(ValueError, "invalid mesh characteristic"):
            MODULE.verify_mesh_family(rows)


if __name__ == "__main__":
    unittest.main()
