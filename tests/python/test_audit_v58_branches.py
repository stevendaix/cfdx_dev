import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from audit_v58_branches import parse_lines


def sample():
    lines = []
    for c in range(3):
        lines.append(f"FPHYS_BRANCH_V58 iteration=1 component={c} faces=634 ambiguous=1 "
                     "effective_kinks=1 unresolved=0 inert=0 flux_crossings=1 roundoff_limited=0 "
                     "sided=0.4 halving=1e-7 certifiable=0 "
                     "native_rc_reclosure_at_offsets=0 physical_qualified=0 committed=0")
        for f in range(3):
            lines.append(f"FPHYS_BRANCH_FACE_V58 iteration=1 component={c} "
                         f"face={f} owner=0 neighbour=1 phi=0 dphi=1 right=2 "
                         "left=5 sided=0.6 halving=0 kink=1 unresolved=0 "
                         "native_rc_reclosure_at_offsets=0 physical_qualified=0 committed=0")
    return lines


class TestV58Proof(unittest.TestCase):
    def test_negative_certificate(self):
        result = parse_lines(sample())
        self.assertEqual(result["effective_kinks"], 3)
        self.assertFalse(result["global_newton_qualified"])

    def test_missing_component(self):
        with self.assertRaises(ValueError): parse_lines(sample()[:4])

    def test_duplicate_group(self):
        with self.assertRaises(ValueError): parse_lines(sample() + sample()[:1])

    def test_forged_qualification(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("physical_qualified=0", "physical_qualified=1") for x in sample()])

    def test_wrong_rc_claim(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("native_rc_reclosure_at_offsets=0",
                                   "native_rc_reclosure_at_offsets=1") for x in sample()])

    def test_unreal_certificate(self):
        with self.assertRaises(ValueError):
            parse_lines([x.replace("certifiable=0", "certifiable=1") for x in sample()])

    def test_truncated_face_detail(self):
        with self.assertRaises(ValueError): parse_lines(sample()[:-1])

    def test_nonfinite(self):
        with self.assertRaises(ValueError): parse_lines([x.replace("sided=0.4", "sided=nan") for x in sample()])


if __name__ == "__main__": unittest.main()
