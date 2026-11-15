"""Fail-closed tests for V76 plots and the V75 budget import.

No synthetic data is ever sent to the reference plotter or reported as
measured CFDX curves. These deliberately fabricated rows test the parser only.
"""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "bfs_comparison_plots_v76.py"
spec = importlib.util.spec_from_file_location("bfsplot", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def row(**kwargs):
    values = dict(iteration=23, retry=8, cycle=2, cell=0, equation="k", trace_selected=0,
        committed=0, convection_delta=0.01, diffusion_delta=-0.04,
        production_delta=0.03, destruction_delta=-0.02,
        cross_delta=0.01, terms_delta_sum=-0.01,
        raw_delta=-0.01, closure_error=0.0)
    values.update(kwargs)
    return "BFS_SST_BUDGET_V75 " + " ".join(f"{k}={v}" for k,v in values.items())


class PlotAndBudgetContract(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as d:
            filename = Path(d) / "synthetic.log"
            filename.write_text(text)
            return module.parse_budgets(filename)

    def test_valid_signed_budget(self):
        records = self.parse(row())
        self.assertEqual(len(records), 1)
        self.assertAlmostEqual(records[0][6], -.01)

    def test_invalid_sum(self):
        with self.assertRaisesRegex(ValueError, "nonconservative"):
            self.parse(row(raw_delta=.02))

    def test_false_commit(self):
        with self.assertRaisesRegex(ValueError, "fabricated commit"):
            self.parse(row(committed=1))

    def test_duplicate(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.parse(row() + "\n" + row())

    def test_nonfinite(self):
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            self.parse(row(production_delta="nan"))

    def test_fake_sum_certificate(self):
        with self.assertRaisesRegex(ValueError, "certificate inconsistent"):
            self.parse(row(terms_delta_sum=123.))

    def test_fake_closure_certificate(self):
        with self.assertRaisesRegex(ValueError, "certificate inconsistent"):
            self.parse(row(closure_error=0.1))

    def test_incoherent_raw_states(self):
        with self.assertRaisesRegex(ValueError, "residual delta inconsistent"):
            self.parse(row(raw_before=1.0, raw_after=2.0))

    def test_no_budget(self):
        with self.assertRaisesRegex(ValueError, "no V75"):
            self.parse("INCOMPRESSIBLE_ITER n=22\n")

    def test_reference_reattachment(self):
        root = module.bzero(module.CF_EXP[:,:2])
        predicted = module.bzero(module.CF_NASA)
        self.assertAlmostEqual(root, 6.278656716417911, places=8)
        self.assertAlmostEqual(predicted, 6.543471646121761, places=8)
        self.assertGreater(predicted, root)

    def test_external_curve_validator(self):
        with tempfile.TemporaryDirectory() as d:
            filename=Path(d)/"curve.csv"
            filename.write_text("x_over_h,cf\n3,-0.001\n5,-0.0007\n7,0.0005\n9,0.001\n")
            values=module.load_code_curve_csv(filename,"cf")
            self.assertEqual(values.shape,(4,2))
            self.assertAlmostEqual(module.bzero(values),6.166666666666667)
            filename.write_text("x_over_h,cf\n3,-0.001\n7,-0.0007\n7,0.0005\n9,0.001\n")
            with self.assertRaisesRegex(ValueError,"strictly"):
                module.load_code_curve_csv(filename,"cf")

    def test_v74_transcription(self):
        self.assertGreater(module.V74["k"][2],module.V74["k"][1])
        self.assertGreater(module.V74["omega"][2],module.V74["omega"][1])
        self.assertLess(module.V74["mass"][2],module.V74["mass"][1])

if __name__ == '__main__':
    unittest.main()