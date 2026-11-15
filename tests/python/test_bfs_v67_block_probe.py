import unittest

from scripts.audit_bfs_v67_block_probe import FRACTIONS, RETRIES, audit


def fake_log():
    parts = []
    for retry in RETRIES:
        beta = 2.0 ** (-retry)
        parts.append(
            f"GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=23 retry={retry:d} "
            f"beta={beta:.12e} ratio={0.99:.12e}"
        )
        for k, w in FRACTIONS:
            ratio = 0.99 if w < 0.75 else 0.998
            growth = 1.1 if k < 1.0 else 2.0
            parts.append(
                f"BFS_BLOCK_COUNTERFACTUAL_V67 iteration=23 retry={retry:d} beta={beta:.12e} "
                f"k_fraction={k:.12e} omega_fraction={w:.12e} merit_ratio={ratio:.12e} "
                f"largest_growth={growth:.12e} worst_block=5 "
                f"joint_ok={int(growth <= 1.2):d} monotone_ok={int(ratio <= 0.995):d} "
                f"original_guards_ok={int(growth <= 1.2 and ratio <= 0.995):d} "
                f"armijo_reference={1.0 - 0.005 * beta:.12e} "
                "counterfactual=1 transport_solved=0 committed=0 qualified=0"
            )
    # Full-six counterfactual sequence from the same retry-zero full trial.
    for power in range(11):
        alpha = 2.0 ** (-power)
        ratio = 0.98 if power == 1 else 0.999
        if power == 0:
            ratio = 0.99  # same as V65 base trial in this synthetic fixture
        growth = 1.1
        parts.append(
            "BFS_FULL_SIX_COUNTERFACTUAL_V67 iteration=23 retry=0 "
            "alpha=%.12e merit_ratio=%.12e largest_growth=%.12e "
            "worst_block=5 joint_ok=1 monotone_ok=%d original_guards_ok=%d "
            "counterfactual=1 carrier_sole_beta=0 committed=0 qualified=0"
            % (alpha, ratio, growth, ratio <= 0.995, ratio <= 0.995)
        )
    return "\n".join(parts)


class ProbeEvidence(unittest.TestCase):
    def test_complete(self):
        evidence = audit(fake_log())
        self.assertEqual(evidence["counterfactual_count"], 36)
        self.assertEqual(evidence["full_six_count"], 11)
        self.assertFalse(evidence["physical_qualification"])
        self.assertGreater(evidence["would_pass_historical_gates"], 0)

    def test_missing(self):
        with self.assertRaises(ValueError):
            audit("\n".join(fake_log().splitlines()[:-1]))

    def test_duplicate(self):
        t = fake_log()
        with self.assertRaises(ValueError):
            audit(t + "\n" + t.splitlines()[1])

    def test_falsified_commit(self):
        with self.assertRaises(ValueError):
            audit(fake_log().replace("committed=0", "committed=1", 1))

    def test_altered_gate(self):
        with self.assertRaises(ValueError):
            audit(fake_log().replace("original_guards_ok=1", "original_guards_ok=0", 1))

    def test_bad_reference(self):
        with self.assertRaises(ValueError):
            audit(
                fake_log().replace(
                    "armijo_reference=9.950000000000e-01",
                    "armijo_reference=9.000000000000e-01",
                    1,
                )
            )


if __name__ == "__main__":
    unittest.main()
