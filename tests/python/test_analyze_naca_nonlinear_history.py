"""Fail-closed regression for accepted NACA nonlinear histories."""
import csv
import importlib.util
import tempfile
from pathlib import Path
import unittest

SRC = Path(__file__).resolve().parents[2] / 'scripts/analyze_naca_nonlinear_history.py'
SPEC = importlib.util.spec_from_file_location('naca_nonlinear_diag', SRC)
assert SPEC and SPEC.loader
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def rows():
    out = []
    for i in range(1, 12):
        row = {k: 0.0 for k in MOD.SCALARS}
        row.update({k: 0 for k in MOD.COUNTS})
        row.update(iteration=i, momentum_rhs_relative=0.03,
                   continuity_normalized=1e-9, mass_local_linf=2e-9,
                   momentum_rhs_x=0.02, momentum_rhs_y=0.9,
                   momentum_balance_x=0.02, momentum_balance_y=0.6 + i*0.01,
                   momentum_worst_cell_y=71, alpha_u=0.6, alpha_p=0.25)
        out.append(row)
    return out


class NonlinearHistoryDiagnosticsTests(unittest.TestCase):
    def test_pressure_converged_is_not_naca_pass(self):
        report = MOD.diagnose(rows())
        self.assertTrue(report['status'].startswith('DIAGNOSTIC_ONLY'))
        self.assertEqual(report['first_pressure_ahead_of_momentum_iteration'], 1)
        self.assertTrue(report['momentum_plateau_last_window'])
        self.assertEqual(report['worst_independent_momentum_balance']['component'], 'y')
        self.assertEqual(report['worst_independent_momentum_balance']['cell'], 71)

    def test_reject_nonfinite_csv_and_nonmonotone_iterations(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'history.csv'
            def save(data):
                with path.open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(data[0]))
                    writer.writeheader()
                    writer.writerows(data)
            bad = rows()
            bad[4]['momentum_rhs_relative'] = float('nan')
            save(bad)
            with self.assertRaisesRegex(ValueError, 'nonfinite'):
                MOD.read_history(path)
            bad = rows()
            bad[2]['iteration'] = 1
            save(bad)
            with self.assertRaisesRegex(ValueError, 'strictly increasing'):
                MOD.read_history(path)

    def test_reject_missing_momentum_component_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'history.csv'
            path.write_text('iteration,momentum_rhs_relative\n1,1e-3\n')
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                MOD.read_history(path)


if __name__ == '__main__':
    unittest.main()
