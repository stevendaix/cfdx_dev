"""Regression tests: fail-closed physical/mesh cross-family comparison."""

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/naca_compare_mesh_families.py"
spec = importlib.util.spec_from_file_location(
    "naca_compare_mesh_families", SCRIPT
)
m = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(m)


def quality(power=1.8, scale=1.0):
    return [
        {
            "n_surface": level,
            "n_radial": 3 * level // 8,
            "radius": 20.0,
            "radial_power": power,
            "minimum_first_layer_spacing": scale / level,
            "mean_first_layer_spacing": 1.1 * scale / level,
            "maximum_cell_edge_aspect_ratio": 15.0 * scale,
            "inverted_or_degenerate_cells": 0,
        }
        for level in (64, 128, 256)
    ]


def physical_report(cd=0.119):
    parts = {}
    for name, fraction in (("cd", 1.0), ("cdp", 0.4), ("cdv", 0.6)):
        parts[name] = {
            "status": "asymptotic_candidate",
            "extrapolated": cd * fraction,
            "GCI_absolute": 0.002 * cd * fraction,
        }
    return {
        "status": "GCI_CANDIDATE",
        "mesh_family": {"status": "COMPATIBLE_WITH_SINGLE_RATIO_GCI"},
        "gci": {"status": "RICHARDSON_CANDIDATE_ONLY", "result": {"quantities": parts}},
        "cases": [{"status": "CONVERGED", "N": n} for n in (64, 128, 256)],
    }


def metadata():
    return {
        "airfoil": "NACA0012", "reynolds": 1000, "angle_deg": 0,
        "chord": 1, "farfield_radius_chords": 20,
        "velocity_boundary": "fixed_farfield", "pressure_boundary": "fixed_farfield",
        "convection_scheme": "TVD_VANLEER", "wall_shear_method": "face_traction",
        "solver_revision": "checked_commit_hash",
    }


class CompareMeshFamilies(unittest.TestCase):
    def test_no_cfd_data_blocks_instead_of_claiming_physics(self):
        report = m.compare(quality(), quality(power=2.0))
        self.assertEqual(report["status"], "BLOCKED")
        self.assertIn("three-grid", report["reasons"][0])

    def test_different_wall_spacing_detected(self):
        diag = m.compare(quality(), quality(power=1.0, scale=13.0))
        self.assertEqual(diag["geometry"]["status"], "DIFFERENT_RESOLUTION")
        self.assertGreater(diag["geometry"]["pairs"][0]["first_wall_layer_ratio"], 12)

    def test_two_converged_three_grid_families_can_agree(self):
        report = m.compare(
            quality(),
            quality(power=2.0),
            physical_report(),
            physical_report(0.1192),
            metadata(),
            metadata(),
        )
        self.assertEqual(report["status"], "CROSS_MESH_GCI_AGREEMENT_CANDIDATE")
        self.assertTrue(report["gci_comparison"]["cdv"]["gci_bands_overlap"])
        self.assertEqual(report["reference"]["physical_validation"], "NOT_PERFORMED")

    def test_total_drag_cancellation_does_not_hide_viscous_error(self):
        a, b = physical_report(), physical_report()
        b["gci"]["result"]["quantities"]["cdv"]["extrapolated"] *= 1.3
        result = m.compare(quality(), quality(), a, b, metadata(), metadata())
        self.assertEqual(result["status"], "BLOCKED")
        self.assertFalse(result["gci_comparison"]["cdv"]["agreement_candidate"])

    def test_nonconverged_family_rejected(self):
        a, b = physical_report(), physical_report()
        b["status"] = "BLOCKED"
        report = m.compare(quality(), quality(), a, b, metadata(), metadata())
        self.assertEqual(report["status"], "BLOCKED")

    def test_unverified_physics_inputs_rejected(self):
        meta = metadata()
        meta["reynolds"] = 2000
        report = m.compare(
            quality(),
            quality(),
            physical_report(),
            physical_report(),
            metadata(),
            meta,
        )
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(any("reynolds" in reason for reason in report["reasons"]))

    def test_incomplete_metadata_rejected(self):
        meta = metadata()
        meta.pop("solver_revision")
        report = m.compare(
            quality(), quality(), physical_report(), physical_report(), meta, metadata()
        )
        self.assertEqual(report["status"], "BLOCKED")

    def test_invalid_values_rejected(self):
        q = quality()
        q[0]["maximum_cell_edge_aspect_ratio"] = float("nan")
        with self.assertRaises(ValueError):
            m.compare(q, quality())

    def test_cli_creates_diagnostic_with_failed_gate(self):
        import subprocess
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            a, b, output = path / "a.json", path / "b.json", path / "report.json"
            a.write_text(json.dumps(quality()), encoding="utf-8")
            b.write_text(json.dumps(quality(scale=13)), encoding="utf-8")
            cmd = [
                sys.executable,
                str(SCRIPT),
                "--a-quality",
                str(a),
                "--b-quality",
                str(b),
                "--output",
                str(output),
            ]
            status = subprocess.run(cmd, capture_output=True, check=False)
            self.assertEqual(status.returncode, 2)
            self.assertEqual(json.loads(output.read_text())["status"], "BLOCKED")


if __name__ == "__main__":
    unittest.main()