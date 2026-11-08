"""Negative controls: never promote wrong operators to validated physical CFD."""

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SOURCE = (
    Path(__file__).resolve().parents[2] / "scripts" / "naca_same_state_tvd_audit.py"
)
SPEC = importlib.util.spec_from_file_location("naca_same_state_tvd_audit", SOURCE)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class SameStateTraceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.face = self.root / "face.csv"
        self.picard = self.root / "picard.csv"
        self.linear = self.root / "linear.csv"
        self.faces = []
        self.picards = []
        self.linears = []
        for it in (1, 2):
            for face, owner, neighbour, mass in ((5, 0, 1, 1.0), (6, 1, -1, 0.2)):
                up = 0.4 if face == 5 else 0.01
                delta = 0.20 if face == 5 else 0.0
                pred = mass + 0.4 if face == 5 else mass
                self.faces.append(
                    {
                        "iteration": it,
                        "face": face,
                        "owner": owner,
                        "neighbour": neighbour,
                        "mass_flux": mass,
                        "upwind_convection": up,
                        "tvd_convection": up + delta,
                        "deferred_tvd_flux": delta,
                        "predictor_hbya_flux": pred,
                        "rc_relaxed_flux": mass,
                        "rc_full_flux": mass,
                        "conservative_flux": mass,
                        "rc_pressure_response": pred - mass,
                        "pressure_relaxation_split": 0,
                        "simplec_extra_flux": 0,
                        "decomposition_error": 0,
                        "internal_pressure_operator": (
                            pred - mass if face == 5 else "nan"
                        ),
                        "pressure_operator_error": 0 if face == 5 else "nan",
                    }
                )
                self.picards.append(
                    {
                        "iteration": it,
                        "face": face,
                        "owner": owner,
                        "neighbour": neighbour,
                        "high_order_fraction": 1,
                        "assembly_mass_flux": mass,
                        "upwind_flux": up,
                        "selected_flux": up + delta,
                        "deferred_tvd_flux": delta,
                    }
                )
            for c in (0, 1, 2):
                self.linears.append(
                    {
                        "iteration": it,
                        "component": c,
                        "alpha_u": 0.6,
                        "linear_tolerance": 1e-10,
                        "relaxed_true_l2": 1e-13,
                        "relaxed_true_relative": 1e-12,
                        "unrelaxed_original_l2": 0.2,
                        "unrelaxed_original_relative": 0.02,
                        "relaxation_displacement_l2": 0.2,
                        "velocity_change_linf": 0.001,
                        "relaxed_true_linf": 1e-13,
                        "original_operator_linf": 0.2,
                        "worst_original_cell": 0,
                    }
                )
        self.write()

    def write(self):
        for path, cols, rows in (
            (self.face, AUDIT.FACE_COLUMNS, self.faces),
            (self.picard, AUDIT.PICARD_COLUMNS, self.picards),
            (self.linear, AUDIT.LINEAR_COLUMNS, self.linears),
        ):
            with path.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=cols)
                writer.writeheader()
                writer.writerows(rows)

    def assess(self):
        return AUDIT.analyze(self.face, self.linear, (0, 1), 1, 2, self.picard)

    def test_valid_operator_trace_is_not_physical_qualification(self):
        d = self.assess()
        self.assertEqual(d["trace_consistency"], "PASS")
        self.assertEqual(d["physical_qualification"], "NOT_EVALUATED")
        self.assertEqual(
            d["max_absolute_closure_errors"]["previous_accepted_flux_lag"], 0
        )

    def test_wrong_deferred_flux_fails(self):
        self.faces[1]["deferred_tvd_flux"] = 0.3
        self.write()
        with self.assertRaisesRegex(ValueError, "deferred"):
            self.assess()

    def test_pressure_operator_mismatch_fails(self):
        self.faces[0]["internal_pressure_operator"] = 0.05
        self.write()
        with self.assertRaisesRegex(ValueError, "independent_internal_pressure"):
            self.assess()

    def test_stale_predictor_mass_flux_fails(self):
        self.picards[2]["assembly_mass_flux"] = 1.1
        self.write()
        with self.assertRaisesRegex(ValueError, "stale/wrong Picard mass flux"):
            self.assess()

    def test_false_linear_certification_fails(self):
        self.linears[3]["relaxed_true_relative"] = 1e-7
        self.write()
        with self.assertRaisesRegex(ValueError, "TRUE residual"):
            self.assess()

    def test_missing_iteration_fails(self):
        self.faces = self.faces[:2]
        self.write()
        with self.assertRaisesRegex(ValueError, "missing iterations"):
            self.assess()

    def test_duplicate_face_fails(self):
        self.faces.append(self.faces[0].copy())
        self.write()
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.assess()

    def test_incomplete_linear_certificate_fails(self):
        self.linears = self.linears[:-1]
        self.write()
        with self.assertRaisesRegex(ValueError, "missing/extra"):
            self.assess()

    def test_invalid_pressure_nan_fails(self):
        self.faces[0]["pressure_operator_error"] = "nan"
        self.write()
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            self.assess()


if __name__ == "__main__":
    unittest.main()
