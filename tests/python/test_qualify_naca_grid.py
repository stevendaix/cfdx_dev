import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "qualify_naca_grid.py"
SPEC = importlib.util.spec_from_file_location("cfdx_qualify_naca_grid", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
QUALIFYMOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QUALIFYMOD)
qualify = QUALIFYMOD.qualify


def row(h):
    cd = 0.12 + 2 * h * h
    return {
        "case": str(h),
        "h": str(h),
        "converged": "yes",
        "res_mom": "1e-10",
        "res_mass": "1e-11",
        "cd": str(cd),
        "cdp": str(0.04 + h * h),
        "cdv": str(0.08 + h * h),
        "cl": "0",
        "cm": "0",
    }


class GridQualificationTests(unittest.TestCase):
    def test_mms_grid_order_two(self):
        q = qualify([row(0.25), row(0.125), row(0.0625)])
        self.assertAlmostEqual(q["quantities"]["cd"]["order"], 2.0)
        self.assertAlmostEqual(q["quantities"]["cd"]["extrapolated"], 0.12)

    def test_unconverged_grid_is_rejected(self):
        rows = [row(0.25), row(0.125), row(0.0625)]
        rows[0]["converged"] = "no"
        with self.assertRaisesRegex(ValueError, "NOT attained"):
            qualify(rows)

    def test_bad_residual_is_rejected(self):
        rows = [row(0.25), row(0.125), row(0.0625)]
        rows[2]["res_mass"] = "1e-5"
        with self.assertRaisesRegex(ValueError, "exceeds"):
            qualify(rows)

    def test_nonmonotone_cd_has_no_gci(self):
        rows = [row(0.25), row(0.125), row(0.0625)]
        rows[2]["cd"] = rows[0]["cd"]
        rows[2]["cdv"] = str(float(rows[2]["cd"]) - float(rows[2]["cdp"]))
        self.assertIsNone(qualify(rows)["quantities"]["cd"]["GCI"])

    def test_bad_force_partition_is_rejected(self):
        rows = [row(0.25), row(0.125), row(0.0625)]
        rows[0]["cdv"] = "0.5"
        with self.assertRaisesRegex(ValueError, "Cd !="):
            qualify(rows)

    def test_non_geometric_grid_family_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unequal refinement"):
            qualify([row(0.25), row(0.11), row(0.0625)])


if __name__ == "__main__":
    unittest.main()
