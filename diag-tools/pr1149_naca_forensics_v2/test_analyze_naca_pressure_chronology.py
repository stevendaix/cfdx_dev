from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "validation"
    / "analyze_naca_pressure_chronology.py"
)
SPEC = importlib.util.spec_from_file_location(
    "analyze_naca_pressure_chronology", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class NacaPressureChronologyTests(unittest.TestCase):
    def test_parse_extracts_solver_and_coupling_fields(self):
        text = (
            "NACA_PRESSURE_CORRECTION iter=1 corrector=0 status=1 "
            "iterations=14 residual=1.2e-11 relative=1.6e-11 rows=12 "
            "nnz=61 rhs_linf=-5.7 rhs_cell=3 correction_linf=5.7 "
            "correction_cell=3 pressure_before_linf=0.0 pressure_cell=0 "
            "diagonal_min=1.0 diagonal_abs_max=2.0 diagonal_cell=1 "
            "minimum_diagonal_cell=2 minimum_diagonal_rhs=0.5 "
            "minimum_diagonal_correction=0.25 rhs_diagonal_ratio_linf=0.75 "
            "rhs_diagonal_ratio_cell=2 diagonal_dominance_min=1.5 "
            "minimum_dominance_cell=4 row_sum_linf=0.0 "
            "row_abs_sum_linf=4.0\n"
        )
        records = MODULE.parse_pressure_correction(text)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record["iter"], 1)
        self.assertEqual(record["status"], 1)
        self.assertEqual(record["iterations"], 14)
        self.assertEqual(record["rhs_cell"], 3)
        self.assertEqual(record["correction_linf"], 5.7)

    def test_first_growth_detects_consecutive_ratio(self):
        lines = []
        for iteration, correction in ((1, 1.0), (2, 1.5), (3, 30.0)):
            lines.append(
                f"NACA_PRESSURE_CORRECTION iter={iteration} corrector=0 "
                f"status=1 iterations=5 residual=1e-12 relative=1e-12 "
                f"rows=4 nnz=12 rhs_linf=1.0 rhs_cell=0 "
                f"correction_linf={correction} correction_cell=0 "
                f"pressure_before_linf=0.0 pressure_cell=0 "
                f"diagonal_min=1.0 diagonal_abs_max=1.0 diagonal_cell=0 "
                f"minimum_diagonal_cell=0 minimum_diagonal_rhs=0.0 "
                f"minimum_diagonal_correction=0.0 "
                f"rhs_diagonal_ratio_linf=1.0 rhs_diagonal_ratio_cell=0 "
                f"diagonal_dominance_min=1.0 minimum_dominance_cell=0 "
                f"row_sum_linf=0.0 row_abs_sum_linf=2.0"
            )
        records = MODULE.parse_pressure_correction("\n".join(lines))
        event = MODULE.first_growth(records, "correction_linf", 10.0)
        self.assertIsNotNone(event)
        self.assertEqual(event["previous_iter"], 2)
        self.assertEqual(event["current_iter"], 3)
        self.assertEqual(event["factor"], 20.0)

    def test_largest_clean_solve_flags_converged_oversized_corrections(self):
        lines = [
            "NACA_PRESSURE_CORRECTION iter=62 corrector=0 status=1 "
            "iterations=14 residual=1.6e-11 relative=1.6e-11 rows=49248 "
            "nnz=1000 rhs_linf=2.97e4 rhs_cell=7 correction_linf=2.97e4 "
            "correction_cell=7 pressure_before_linf=1.0 pressure_cell=7 "
            "diagonal_min=1.0 diagonal_abs_max=2.0 diagonal_cell=1 "
            "minimum_diagonal_cell=2 minimum_diagonal_rhs=0.5 "
            "minimum_diagonal_correction=0.25 rhs_diagonal_ratio_linf=9.9 "
            "rhs_diagonal_ratio_cell=7 diagonal_dominance_min=1.5 "
            "minimum_dominance_cell=4 row_sum_linf=0.0 row_abs_sum_linf=4.0",
            "NACA_PRESSURE_CORRECTION iter=63 corrector=0 status=2 "
            "iterations=1000 residual=1e300 relative=1e300 rows=49248 "
            "nnz=1000 rhs_linf=1e300 rhs_cell=7 correction_linf=1e300 "
            "correction_cell=7 pressure_before_linf=1e300 pressure_cell=7 "
            "diagonal_min=1.0 diagonal_abs_max=2.0 diagonal_cell=1 "
            "minimum_diagonal_cell=2 minimum_diagonal_rhs=0.5 "
            "minimum_diagonal_correction=0.25 rhs_diagonal_ratio_linf=9.9 "
            "rhs_diagonal_ratio_cell=7 diagonal_dominance_min=1.5 "
            "minimum_dominance_cell=4 row_sum_linf=0.0 row_abs_sum_linf=4.0",
        ]
        records = MODULE.parse_pressure_correction("\n".join(lines))
        clean = MODULE.largest_clean_solve(records)
        self.assertIsNotNone(clean)
        self.assertEqual(clean["iter"], 62)
        self.assertEqual(clean["correction_linf"], 2.97e4)
        self.assertEqual(clean["relative"], 1.6e-11)

    def test_cli_reports_records_and_verdict(self):
        line = (
            "NACA_PRESSURE_CORRECTION iter=1 corrector=0 status=1 "
            "iterations=3 residual=1e-12 relative=1e-12 rows=4 nnz=12 "
            "rhs_linf=1.0 rhs_cell=0 correction_linf=1.0 correction_cell=0 "
            "pressure_before_linf=0.0 pressure_cell=0 diagonal_min=1.0 "
            "diagonal_abs_max=1.0 diagonal_cell=0 minimum_diagonal_cell=0 "
            "minimum_diagonal_rhs=0.0 minimum_diagonal_correction=0.0 "
            "rhs_diagonal_ratio_linf=1.0 rhs_diagonal_ratio_cell=0 "
            "diagonal_dominance_min=1.0 minimum_dominance_cell=0 "
            "row_sum_linf=0.0 row_abs_sum_linf=2.0\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "chronology.log"
            log.write_text(line, encoding="utf-8")
            import contextlib
            import io
            import sys

            captured = io.StringIO()
            with contextlib.redirect_stdout(captured):
                previous_argv = sys.argv
                sys.argv = [MODULE.__name__, str(log)]
                try:
                    exit_code = MODULE.main()
                finally:
                    sys.argv = previous_argv
        self.assertEqual(exit_code, 0)
        self.assertIn("records=1", captured.getvalue())
        self.assertIn("LARGEST_CLEAN_SOLVE", captured.getvalue())


if __name__ == "__main__":
    unittest.main()
