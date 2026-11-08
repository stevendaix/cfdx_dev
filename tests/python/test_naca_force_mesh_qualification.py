"""Independent negative and positive contracts for NACA force+mesh evidence."""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import naca_force_mesh_qualification as m
from naca_surface_force_oracle import COLUMNS


def mesh(similar=True):
    return [
        {
            "n_surface": n,
            "n_radial": n // 2,
            "radius": 20.0,
            "radial_power": 1.8,
            "minimum_first_layer_spacing": (1 / n if similar else 1 / (n**1.6)),
            "inverted_or_degenerate_cells": 0,
            "radial_refinement_strategy": "nested_base_24",
            "mean_first_layer_spacing": (1.25 / n if similar else 1.25 / n**1.6),
            "maximum_first_layer_spacing": (1.5 / n if similar else 1.5 / n**1.6),
            "minimum_absolute_cell_area": 0.001 / n,
            "maximum_cell_edge_aspect_ratio": 8.0,
            "minimum_radial_growth_ratio": 1.0,
            "maximum_radial_growth_ratio": 1.2,
            "maximum_nonorthogonality_degrees": 45.0,
            "p95_nonorthogonality_degrees": 20.0,
        }
        for n in m.LEVELS
    ]


class NacaForceMeshTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.folder = Path(self.tmp.name)
        self.certs = {}
        self.histories = {}
        self.raws = {}
        for n, cd in zip(m.LEVELS, [0.11, 0.13, 0.138]):
            c = self.folder / f"cert_{n}.csv"
            with c.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=m.REQUIRED)
                writer.writeheader()
                writer.writerow(
                    {
                        "case": f"N{n}",
                        "h": 1 / n,
                        "converged": "yes",
                        "res_mom": 1e-11,
                        "res_mass": 2e-11,
                        "cd": cd,
                        "cdp": cd * 0.4,
                        "cdv": cd * 0.6,
                        "cl": 0,
                        "cm": 0,
                    }
                )
            h = self.folder / f"history_{n}.csv"
            h.write_text(
                "iteration,momentum_rhs_relative,continuity_normalized,"
                "mass_local_linf,nonfinite_mass_faces,nonfinite_velocities,"
                "high_order_fraction,face_momentum_x_relative,"
                "face_momentum_y_relative,face_momentum_z_relative,"
                "face_momentum_x_linf,face_momentum_y_linf,face_momentum_z_linf,"
                "mass_normalized_imbalance,convection_scheme,"
                "momentum_x_rhs_relative,"
                "momentum_y_rhs_relative,momentum_z_rhs_relative,"
                "momentum_component_max,pressure_residual,velocity_change_inf,"
                "pressure_change_inf\n"
                "1,1e-11,1e-12,2e-11,0,0,1,1e-11,1e-11,0,1e-11,1e-11,0,"
                "1e-13,TVD,1e-11,1e-11,0,1e-11,1e-12,1e-12,1e-12\n",
                encoding="utf-8",
            )
            self.certs[n] = c
            self.histories[n] = h
            # Synthetic closed 2-D face contour, only for certificate checks.
            # x-normal pressure p=(0.4*Cd)*x; uniform synthetic tau_n_x
            # contributes 0.6*Cd over total area 4; origin y=0 yields Cm=0.
            raw = self.folder / f"raw_{n}.csv"
            with raw.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=COLUMNS)
                writer.writeheader()
                face = 0
                for wall in ("left", "right", "bottom", "top"):
                    for i in range(8):
                        t = (i + 0.5) / 8
                        x, y, nx, ny = {
                            "left": (0, t - 0.5, -1, 0),
                            "right": (1, t - 0.5, 1, 0),
                            "bottom": (t, -0.5, 0, -1),
                            "top": (t, 0.5, 0, 1),
                        }[wall]
                        writer.writerow(
                            {
                                "face": face,
                                "x": x,
                                "y": y,
                                "z": 0.5,
                                "nx": nx,
                                "ny": ny,
                                "nz": 0,
                                "area": 0.125,
                                "pressure": cd * 0.4 * x,
                                "tau_n_x": -cd * 0.6 / 4,
                                "tau_n_y": 0,
                                "tau_n_z": 0,
                                "dynamic_pressure": 1,
                                "reference_area": 1,
                                "reference_chord": 1,
                                "x_ref": 0.25,
                                "y_ref": 0,
                            }
                        )
                        face += 1
            self.raws[n] = raw

    def test_transverse_component_blocks_false_global_convergence(self):
        path = self.histories[128]
        contents = path.read_text()
        path.write_text(
            contents.replace(",TVD,1e-11,1e-11,0,1e-11,", ",TVD,1e-11,2e-5,0,2e-5,")
        )
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any("transverse momentum" in reason for reason in report["reasons"])
        )

    def test_missing_iteration_blocks_incomplete_history(self):
        path = self.histories[128]
        data = path.read_text()
        path.write_text(data.replace("\n1,", "\n2,"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any(
                "incomplete or nonconsecutive" in reason for reason in report["reasons"]
            )
        )

    def test_valid_evidence_produces_only_richardson_candidate(self):
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "GCI_CANDIDATE")
        self.assertEqual(report["gci"]["validation_reference"], "NOT_PERFORMED")

    def test_force_certificate_mass_is_maximum_of_three_measures(self):
        # C++ writes max(continuity, global imbalance, local L_inf).
        # A valid solution must not be rejected when continuity dominates.
        path = self.histories[64]
        with path.open(newline="") as stream:
            row = next(csv.DictReader(stream))
        row["continuity_normalized"] = "9e-9"
        row["mass_normalized_imbalance"] = "7e-9"
        row["mass_local_linf"] = "2e-9"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        cert = self.certs[64]
        with cert.open(newline="") as stream:
            cert_row = next(csv.DictReader(stream))
        cert_row["res_mass"] = "9e-9"
        with cert.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=m.REQUIRED)
            writer.writeheader()
            writer.writerow(cert_row)
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "GCI_CANDIDATE")

    def test_mixed_high_order_schemes_cannot_claim_gci(self):
        # All three meshes would individually pass the residual certificate,
        # but a Richardson ratio across different operators is invalid.
        path = self.histories[128]
        path.write_text(path.read_text().replace(",TVD,", ",TVD_VANLEER,"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any(
                "mixes different convection schemes" in msg for msg in report["reasons"]
            )
        )

    def test_muscl_candidate_remains_experimental_and_needs_three_grids(self):
        for level in m.LEVELS:
            path = self.histories[level]
            path.write_text(path.read_text().replace(",TVD,", ",MUSCL_CELL_LIMITED,"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "GCI_CANDIDATE")
        self.assertEqual(report["convection_scheme"], "MUSCL_CELL_LIMITED")
        self.assertTrue(report["experimental_reconstruction"])
        self.histories[256].unlink()
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")

    def test_mixed_scheme_within_one_history_blocks(self):
        path = self.histories[64]
        path.write_text(
            path.read_text() + "2,1e-11,1e-12,2e-11,0,0,1,1e-11,1e-11,0,1e-11,1e-11,0,"
            "1e-13,TVD_VANLEER,1e-11,1e-11,0,1e-11,1e-12,1e-12,1e-12\n"
        )
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any("mixes convection schemes" in msg for msg in report["reasons"])
        )

    def test_early_upwind_then_final_tvd_is_blocked_by_scheme_consistency(self):
        # Superseded contract (2026-11-08 lot evidence): a nonlinear history
        # that mixes convection schemes is refused, even when only the final
        # accepted operator is used for the certificate.
        path = self.histories[64]
        lines = path.read_text().splitlines()
        earlier = lines[1].replace("1,1e-11,", "1,1e-4,", 1)
        earlier = earlier.replace(",TVD,", ",UPWIND,")
        final = "2" + lines[1][1:]
        path.write_text("\n".join((lines[0], earlier, final)) + "\n")
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any("mixes convection schemes" in msg for msg in report["reasons"])
        )

    def test_missing_raw_surface_cannot_promote_solver_forces(self):
        self.raws[64].unlink()
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(any("raw_64" in reason for reason in report["reasons"]))

    def test_false_force_coefficient_is_detected_by_raw_area_integration(self):
        path = self.raws[64]
        data = path.read_text()
        # Change the pressure on a positive-x face while preserving the
        # certificate and residual history: only raw integration sees it.
        path.write_text(data.replace("0.044000000000000004", "0.084"))
        # Use an unambiguous numeric mutation even if float formatting varies.
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        rows[8]["pressure"] = "3.0"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")

    def test_incompatible_near_wall_refinement_refuses_gci(self):
        report = m.certify(mesh(False), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertEqual(report["gci"]["status"], "NOT_QUALIFIED")

    def test_missing_n256_certificate_cannot_be_promoted(self):
        self.certs[256].unlink()
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertEqual(report["cases"][2]["status"], "NOT_QUALIFIED")

    def test_stale_incorrect_certificate_fails(self):
        self.certs[128].write_text(self.certs[128].read_text().replace("N128", "N129"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(any("wrong mesh label" in item for item in report["reasons"]))

    def test_mismatched_physical_history_fails(self):
        path = self.histories[64]
        original = path.read_text()
        path.write_text(
            original.replace("1,1e-11,1e-12", "1,8e-9,1e-12").replace(
                ",TVD,1e-11,1e-11,0,1e-11,", ",TVD,8e-9,1e-11,0,8e-9,"
            )
        )
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(any("stale" in item for item in report["reasons"]))

    def test_non_monotone_drag_is_not_richardson(self):
        path = self.certs[256]
        path.write_text(path.read_text().replace("0.138", "0.115"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")

    def test_total_drag_gci_cannot_hide_nonmonotone_pressure_and_viscous_components(
        self,
    ):
        # Cd stays on an exact quadratic trend while Cd_pressure reverses.
        # A valid total Cd alone must not certify a pressure+skin-friction GCI.
        path = self.certs[128]
        with path.open(newline="") as stream:
            record = next(iter(csv.DictReader(stream)))
        record["cdp"] = "0.06"
        record["cdv"] = str(float(record["cd"]) - float(record["cdp"]))
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=m.REQUIRED)
            writer.writeheader()
            writer.writerow(record)
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(any("cdp:" in reason for reason in report["reasons"]))

    def test_independent_face_residual_blocks_false_matrix_convergence(self):
        path = self.histories[128]
        old = path.read_text()
        path.write_text(
            old.replace(
                "1,1e-11,1e-12,2e-11,0,0,1,1e-11", "1,1e-11,1e-12,2e-11,0,0,1,0.002"
            )
        )
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any(
                "independent local face momentum" in reason
                for reason in report["reasons"]
            )
        )

    def test_missing_independent_residual_schema_blocks_gci(self):
        path = self.histories[64]
        path.write_text(
            "iteration,momentum_rhs_relative,continuity_normalized,"
            "mass_local_linf,nonfinite_mass_faces,nonfinite_velocities\n"
            "1,1e-11,1e-12,2e-11,0,0\n"
        )
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any("incomplete physical history" in reason for reason in report["reasons"])
        )

    def test_first_order_upwind_cannot_claim_final_high_order_operator(self):
        path = self.histories[64]
        path.write_text(path.read_text().replace(",TVD,", ",UPWIND,"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")
        self.assertTrue(
            any("bounded second-order TVD" in reason for reason in report["reasons"])
        )

    def test_incomplete_final_convection_operator_blocks_gci(self):
        path = self.histories[256]
        path.write_text(
            path.read_text().replace(
                "1,1e-11,1e-12,2e-11,0,0,1,1e-11", "1,1e-11,1e-12,2e-11,0,0,0.5,1e-11"
            )
        )
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")

    def test_negative_certificate_residual_is_not_a_small_residual(self):
        path = self.certs[64]
        path.write_text(path.read_text().replace("1e-11", "-1e-11"))
        history = self.histories[64]
        history.write_text(history.read_text().replace("1e-11", "-1e-11"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")

    def test_negative_continuity_in_history_cannot_be_certified(self):
        path = self.histories[128]
        path.write_text(path.read_text().replace("1e-12", "-1e-12"))
        report = m.certify(mesh(), self.certs, self.histories, self.raws)
        self.assertEqual(report["status"], "BLOCKED")


if __name__ == "__main__":
    unittest.main()
