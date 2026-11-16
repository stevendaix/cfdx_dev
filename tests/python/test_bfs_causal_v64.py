import unittest

from scripts.analyze_bfs_causal_v64 import parse

SUMMARY = (
    "BFS_CAUSAL_V64 iteration=1 retry=0 cells=152 flow_changed=1 "
    "sst_changed=1 relative_flow_delta=1e-2 relative_sst_delta=2e-2 "
    "recomposition_error=1e-12 rc00=1e-15 rc10=1e-15 "
    "rc01=1e-15 rc11=1e-15 qualified=0 committed=0 "
    "four_corners_reassembled=1 phi_flow=0.1 phi_sst=0.01 "
    "phi_interaction=0.02 phi_interaction_face=5 "
    "nut_flow=0.001 nut_sst=0.003 nut_interaction=0.002 "
    "nut_interaction_cell=1"
)


def lines():
    return "\n".join(
        [SUMMARY]
        + [
            f"BFS_CAUSAL_BLOCK_V64 iteration=1 block={i} raw00=1 "
            f"raw10=0.8 raw01=1.1 raw11=0.9 flow_only_delta=0.2 "
            f"sst_only_delta=0.1 interaction=0.05 interaction_relative=0.045 "
            f"worst_after_cell={i} worst_interaction_cell={i} "
            f"qualified=0 committed=0"
            for i in range(6)
        ]
    )


class TestBfsCausalV64(unittest.TestCase):
    def test_six_valid(self):
        d = parse(lines())
        self.assertEqual(len(d["iterations"][0]["blocks"]), 6)
        self.assertFalse(d["physical_qualification"])

    def test_missing(self):
        with self.assertRaises(ValueError):
            parse("")

    def test_truncated_block(self):
        with self.assertRaises(ValueError):
            parse(lines().rsplit("\n", 1)[0])

    def test_forged_commit(self):
        with self.assertRaises(ValueError):
            parse(lines().replace("qualified=0", "qualified=1"))

    def test_unclosed_rc(self):
        with self.assertRaises(ValueError):
            parse(lines().replace("rc11=1e-15", "rc11=1e-5"))

    def test_nonfinite(self):
        with self.assertRaises(ValueError):
            parse(lines().replace("interaction=0.05", "interaction=nan"))

    def test_wrong_cell(self):
        with self.assertRaises(ValueError):
            parse(lines().replace("worst_after_cell=0", "worst_after_cell=152"))

    def test_spurious_block(self):
        with self.assertRaises(ValueError):
            parse(lines() + "\n" + lines().splitlines()[-1])


if __name__ == "__main__":
    unittest.main()
