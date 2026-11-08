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
    return {
        "n_surface": surface,
        "n_radial": radial,
        "radius": 20.0,
        "radial_power": 1.8,
        "minimum_first_layer_spacing": spacing,
        "inverted_or_degenerate_cells": 0,
        "radial_refinement_strategy": "nested_base_24",
        "mean_first_layer_spacing": spacing * 1.25,
        "maximum_first_layer_spacing": spacing * 1.5,
        "minimum_absolute_cell_area": spacing * 0.001,
        "maximum_cell_edge_aspect_ratio": 8.0,
        "minimum_radial_growth_ratio": 1.0,
        "maximum_radial_growth_ratio": 1.2,
        "maximum_nonorthogonality_degrees": 45.0,
        "p95_nonorthogonality_degrees": 20.0,
    }


class MeshFamilyTests(unittest.TestCase):
    def test_geometrically_similar_mesh_family_is_allowed(self):
        self.assertEqual(
            MODULE.verify_mesh_family(
                [row(64, 24, 0.04), row(128, 48, 0.02), row(256, 96, 0.01)]
            )["status"],
            "COMPATIBLE_WITH_SINGLE_RATIO_GCI",
        )

    def test_different_te_sampling_blocks_a_mixed_mesh_family(self):
        rows = [row(64, 24, 0.04), row(128, 48, 0.02), row(256, 96, 0.01)]
        rows[0]["trailing_edge_regularization"] = 0.04
        rows[1]["trailing_edge_regularization"] = 0.04
        rows[2]["trailing_edge_regularization"] = 0.0
        result = MODULE.verify_mesh_family(rows)
        self.assertEqual(result["status"], "BLOCKED")
        self.assertIn("different_trailing_edge_sampling", result["reasons"])

    def test_wall_normal_refinement_mismatch_blocks_gci(self):
        report = MODULE.verify_mesh_family(
            [row(64, 24, 0.05), row(128, 48, 0.0155), row(256, 96, 0.0047)]
        )
        self.assertEqual(report["status"], "BLOCKED")
        self.assertIn("first_layer_distribution_not_similar", report["reasons"])

    def test_bad_geometry_and_inverted_cells_block(self):
        rows = [row(64, 24, 0.04), row(128, 48, 0.02), row(256, 96, 0.01)]
        rows[1]["radius"] = 10.0
        rows[2]["inverted_or_degenerate_cells"] = 1
        self.assertEqual(MODULE.verify_mesh_family(rows)["status"], "BLOCKED")

    def test_mixed_radial_distributions_are_not_comparable(self):
        rows = [row(64, 24, 0.04), row(128, 48, 0.02), row(256, 96, 0.01)]
        rows[0]["radial_refinement_strategy"] = "powerlaw_resampled"
        rows[1]["radial_refinement_strategy"] = "nested_base_24"
        rows[2]["radial_refinement_strategy"] = "nested_base_24"
        result = MODULE.verify_mesh_family(rows)
        self.assertEqual(result["status"], "BLOCKED")
        self.assertIn(
            "inconsistent_or_missing_radial_refinement_strategy", result["reasons"]
        )

    def test_worsening_nonorthogonality_blocks_false_single_ratio_gci(self):
        rows = [row(64, 24, 0.04), row(128, 48, 0.02), row(256, 96, 0.01)]
        rows[1]["maximum_nonorthogonality_degrees"] = 61.0
        rows[2]["maximum_nonorthogonality_degrees"] = 73.0
        result = MODULE.verify_mesh_family(rows)
        self.assertEqual(result["status"], "BLOCKED")
        self.assertIn(
            "nonorthogonality_distribution_changes_with_refinement",
            result["reasons"],
        )

    def test_wall_mean_distribution_fails_even_when_minimum_matches(self):
        rows = [row(64, 24, 0.04), row(128, 48, 0.02), row(256, 96, 0.01)]
        rows[1]["mean_first_layer_spacing"] *= 2
        result = MODULE.verify_mesh_family(rows)
        self.assertIn("first_layer_distribution_not_similar", result["reasons"])

    def test_missing_distribution_metrics_fail_closed(self):
        rows = [row(64, 24, 0.04), row(128, 48, 0.02), row(256, 96, 0.01)]
        del rows[2]["p95_nonorthogonality_degrees"]
        with self.assertRaises(KeyError):
            MODULE.verify_mesh_family(rows)

    def test_duplicate_levels_fail_closed(self):
        rows = [row(64, 24, 0.04), row(64, 24, 0.04), row(256, 96, 0.01)]
        self.assertEqual(MODULE.verify_mesh_family(rows)["status"], "BLOCKED")

    def test_nonfinite_spacing_fails_closed(self):
        rows = [row(64, 24, 0.04), row(128, 48, float("nan")), row(256, 96, 0.01)]
        with self.assertRaisesRegex(ValueError, "invalid mesh characteristic"):
            MODULE.verify_mesh_family(rows)


if __name__ == "__main__":
    unittest.main()
