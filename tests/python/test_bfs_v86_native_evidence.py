"""Fail-closed tests for native V85 evidence inspection and V86 qualification flags."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.audit_bfs_v86_native_100 import read_evidence

A = "GLOBAL_PROGRESS_TRIAL_V85 iteration=1 accepted=1 reason=baseline qualified=0"
B = "GLOBAL_PROGRESS_COMMITTED_V85 iteration=1 committed=1 U=0.99 U_global=0.58 p_continuity=1e-14 k=0.28 omega=0.003 merit=0.99 k_debt=0 omega_debt=0 converged=0 qualified=0"
C = "BFS_FIELD_COMMITTED_V85 iteration=1 cells=20540 reverse_ux_cells=0 ux_min=36 ux_max=50 pressure_min_pa=-10000 pressure_max_pa=-100 pressure_span_pa=9900 converged=0_or_pending"


class TestNativeEvidence(unittest.TestCase):
    def test_accepted_finite(self):
        c, f, t = read_evidence("\n".join((A, B, C)))
        self.assertEqual(list(c), [1])
        self.assertEqual(int(f[1]["reverse_ux_cells"]), 0)

    def test_missing_field_denied(self):
        with self.assertRaisesRegex(ValueError, "matching field"):
            read_evidence("\n".join((A, B)))

    def test_false_certification_denied(self):
        with self.assertRaisesRegex(ValueError, "falsely claims"):
            read_evidence("\n".join((A, B.replace("qualified=0", "qualified=1"), C)))

    def test_negative_velocity_count(self):
        with self.assertRaisesRegex(ValueError, "count"):
            read_evidence(
                "\n".join(
                    (A, B, C.replace("reverse_ux_cells=0", "reverse_ux_cells=20541"))
                )
            )

    def test_forged_pressure_range(self):
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            read_evidence(
                "\n".join(
                    (A, B, C.replace("pressure_span_pa=9900", "pressure_span_pa=9000"))
                )
            )

    def test_no_accepted_trial(self):
        with self.assertRaisesRegex(ValueError, "lacks an accepted"):
            read_evidence("\n".join((A.replace("accepted=1", "accepted=0"), B, C)))

    def test_duplicate_committed(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            read_evidence("\n".join((A, B, B, C)))

    def test_nonfinite_refused(self):
        with self.assertRaisesRegex(ValueError, "Bad native physical"):
            read_evidence("\n".join((A, B.replace("k=0.28", "k=nan"), C)))

    def test_no_fabricated_states(self):
        c, f, t = read_evidence(
            "BFS_PHASE_V10 iteration=1 corrector=0 phase=PRESSURE_LINEAR_BEGIN"
        )
        self.assertFalse(c)
        self.assertFalse(f)

    def test_bad_min_velocity_cell(self):
        with self.assertRaisesRegex(ValueError, "worst-velocity"):
            read_evidence("\n".join((A, B, C + " ux_min_cell=20540 speed_max=62.0")))

    def test_bad_speed_max(self):
        with self.assertRaisesRegex(ValueError, "maximum speed"):
            read_evidence("\n".join((A, B, C + " ux_min_cell=3 speed_max=nan")))


if __name__ == "__main__":
    unittest.main()
