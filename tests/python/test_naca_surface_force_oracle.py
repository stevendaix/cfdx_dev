"""Exact synthetic surface integrals and deliberate corruptions, no CFD solve."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import naca_surface_force_oracle as oracle


def square_faces(pressure_offset=0.0):
    """Unit-square contour: p=x, tau*n=(-0.1,0,0), F=(1.4,0), Mz=-0.7."""
    result = []
    index = 0
    for wall in ("left", "right", "bottom", "top"):
        for k in range(8):
            t = (k + 0.5) / 8
            x, y, nx, ny = {
                "left": (0, t, -1, 0),
                "right": (1, t, 1, 0),
                "bottom": (t, 0, 0, -1),
                "top": (t, 1, 0, 1),
            }[wall]
            result.append(
                {
                    "face": index,
                    "x": x,
                    "y": y,
                    "z": 0.5,
                    "nx": nx,
                    "ny": ny,
                    "nz": 0,
                    "area": 0.125,
                    "pressure": x + pressure_offset,
                    "tau_n_x": -0.1,
                    "tau_n_y": 0,
                    "tau_n_z": 0,
                    "dynamic_pressure": 1,
                    "reference_area": 1,
                    "reference_chord": 1,
                    "x_ref": 0.25,
                    "y_ref": 0,
                }
            )
            index += 1
    return result


class NacaSurfaceOracleTests(unittest.TestCase):
    def test_exact_pressure_drag_shear_drag_and_pitch_moment(self):
        result = oracle.integrate(square_faces())
        self.assertEqual(result["status"], "INDEPENDENT_SURFACE_INTEGRAL_ONLY")
        self.assertAlmostEqual(result["coefficients"]["cdp"], 1.0, places=12)
        self.assertAlmostEqual(result["coefficients"]["cdv"], 0.4, places=12)
        self.assertAlmostEqual(result["coefficients"]["cd"], 1.4, places=12)
        self.assertAlmostEqual(result["coefficients"]["cl"], 0, places=12)
        self.assertAlmostEqual(result["coefficients"]["cm"], -0.7, places=12)

    def test_uniform_pressure_gauge_cannot_change_force(self):
        before = oracle.integrate(square_faces())["coefficients"]
        after = oracle.integrate(square_faces(pressure_offset=52.7))["coefficients"]
        for key in before:
            self.assertAlmostEqual(before[key], after[key], places=11)

    def test_closed_face_sum_agrees_with_certificate(self):
        cert = {"cd": "1.4", "cdp": "1", "cdv": ".4", "cl": "0", "cm": "-.7"}
        self.assertEqual(
            oracle.integrate(square_faces(), cert)["status"],
            "INDEPENDENT_SURFACE_FORCE_CONSISTENT",
        )

    def test_wrong_certificate_drag_is_detected(self):
        with self.assertRaisesRegex(ValueError, "cd: independent surface integral"):
            oracle.integrate(
                square_faces(),
                {"cd": "0.9", "cdp": "1", "cdv": ".4", "cl": "0", "cm": "-.7"},
            )

    def test_missing_wall_face_detects_geometric_nonclosure(self):
        with self.assertRaisesRegex(ValueError, "incomplete or incorrectly oriented"):
            oracle.integrate(square_faces()[:-1])

    def test_reversed_pressure_normal_detects_bad_orientation(self):
        rows = square_faces()
        for row in rows[:8]:
            row["nx"] *= -1
        with self.assertRaisesRegex(ValueError, "incomplete or incorrectly oriented"):
            oracle.integrate(rows)

    def test_changing_reference_area_on_one_face_is_rejected(self):
        rows = square_faces()
        rows[0]["reference_area"] = 2
        with self.assertRaisesRegex(ValueError, "inconsistent reference normalization"):
            oracle.integrate(rows)

    def test_nonfinite_face_data_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / "raw.csv"
            rows = square_faces()
            rows[0]["pressure"] = float("nan")
            self.write_csv(path, rows)
            with self.assertRaisesRegex(ValueError, "nonfinite face data"):
                oracle.read_rows(path)

    def test_duplicate_face_id_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / "raw.csv"
            rows = square_faces()
            rows[1]["face"] = rows[0]["face"]
            self.write_csv(path, rows)
            with self.assertRaisesRegex(ValueError, "duplicate wall face"):
                oracle.read_rows(path)

    def test_nonunit_normal_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / "raw.csv"
            rows = square_faces()
            rows[0]["nx"] = 0.1
            self.write_csv(path, rows)
            with self.assertRaisesRegex(ValueError, "unit vector"):
                oracle.read_rows(path)

    def test_failed_cli_clears_stale_success(self):
        with tempfile.TemporaryDirectory() as t:
            directory = Path(t)
            raw = directory / "raw.csv"
            output = directory / "report.json"
            self.write_csv(raw, square_faces())
            output.write_text('{"status":"PASS"}\n')
            raw.write_text("bad,header\n1,2\n")
            cmd = [
                sys.executable,
                str(ROOT / "scripts/naca_surface_force_oracle.py"),
                "--raw-faces",
                str(raw),
                "--output",
                str(output),
            ]
            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertFalse(output.exists())

    @staticmethod
    def write_csv(path, rows):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=oracle.COLUMNS)
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
