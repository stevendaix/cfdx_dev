from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "validation" / "analyze_naca_trace.py"
SPEC = importlib.util.spec_from_file_location("analyze_naca_trace", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class NacaTraceAnalyzerTests(unittest.TestCase):
    def test_parse_trace_record(self):
        records = MODULE.parse_trace(
            "prefix NACA_MOMENTUM_TRACE seq=0 iter=1 max_p=2.5 max_gradp=3 "
            "max_U=1 max_phi=4 max_divphi=0.1 max_mu_eff=0.001 "
            "boundary_flux_sum=0 boundary_flux_abs_sum=1 nonfinite_total=0"
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["seq"], 0)
        self.assertAlmostEqual(records[0]["max_gradp"], 3.0)

    def test_detect_first_consecutive_growth(self):
        records = [
            {"seq": 0, "iter": 1, "max_p": 1.0},
            {"seq": 1, "iter": 1, "max_p": 12.0},
        ]
        event = MODULE.first_growth(records, 10.0)
        self.assertIsNotNone(event)
        self.assertEqual(event["metric"], "max_p")
        self.assertEqual(event["current_seq"], 1)

    def test_no_growth_is_not_a_failure(self):
        records = [
            {"seq": 0, "iter": 1, "max_U": 1.0},
            {"seq": 1, "iter": 2, "max_U": 2.0},
        ]
        self.assertIsNone(MODULE.first_growth(records, 10.0))

    def test_nonfinite_values_are_ignored_by_ratio_detector(self):
        records = [
            {"seq": 0, "iter": 1, "max_U": float("nan")},
            {"seq": 1, "iter": 2, "max_U": 100.0},
        ]
        self.assertIsNone(MODULE.first_growth(records, 10.0))


if __name__ == "__main__":
    unittest.main()
