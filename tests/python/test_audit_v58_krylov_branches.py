"""Mutation tests: the real V58 Krylov failure must remain fail-closed."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from audit_v58_krylov_branches import parse_lines


FIXTURE = Path(__file__).resolve().parents[1] / "data" / "bfs152_v58_krylov_evidence.log"


def original():
    return FIXTURE.read_text().splitlines()


class TestV58KrylovFailures(unittest.TestCase):
    def test_real_bfs152_blocked(self):
        record = parse_lines(original())
        self.assertEqual(record["status"], "BLOCKED_KRYLOV_NONADDITIVE")
        self.assertEqual(record["face_ids_root"], [5, 65])
        self.assertEqual(record["near_triple_tie_face_ids"], [10, 593])
        self.assertFalse(record["global_newton_qualified"])

    def test_missing_krylov(self):
        with self.assertRaises(ValueError):
            parse_lines([x for x in original() if "direction=6 " not in x])

    def test_missing_branch_sample(self):
        with self.assertRaises(ValueError):
            parse_lines([x for x in original() if not (
                x.startswith("FPHYS_KRYLOV_BRANCH_V58 ") and "sample=17 " in x)])

    def test_forged_success(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("qualified=0", "qualified=1") if
                         x.startswith("FPHYS_JFNK_V52 ") else x for x in original()])

    def test_forged_rhie_chow(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("native_rc_reclosure_at_offsets=0",
                                   "native_rc_reclosure_at_offsets=1") for x in original()])

    def test_forged_commit(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("production_committed=0", "production_committed=1")
                         for x in original()])

    def test_missing_face_metadata(self):
        with self.assertRaises(ValueError):
            parse_lines([x for x in original() if not
                         (x.startswith("FPHYS_LIMITER_STATE_V58 ") and "face=5 " in x)])

    def test_corrupt_krylov_residual(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("linear_true_relative=9.721417820905e-01",
                                   "linear_true_relative=1.000000000000e-02")
                         for x in original()])

    def test_nonfinite_face(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("phi=4.277410800000e-03", "phi=nan", 1)
                         if x.startswith("FPHYS_KRYLOV_FACE_V58 ") else x
                         for x in original()])

    def test_duplicate_jfnk_result(self):
        with self.assertRaises(ValueError):
            lines = original()
            parse_lines(lines + [next(x for x in lines if x.startswith("FPHYS_JFNK_V52 "))])


if __name__ == "__main__":
    unittest.main()
