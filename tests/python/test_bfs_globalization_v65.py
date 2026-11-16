import unittest

from scripts.audit_bfs_globalization_v65 import audit


class V65AuditTests(unittest.TestCase):
    @staticmethod
    def fixture():
        t = "GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=1 retry=0 beta=1 baseline=1 old_merit=3 trial_merit=2 ratio=0.666 largest_growth=1 worst_block=0 joint_ok=1 monotone_ok=1 accepted=1 old_rc_gap=1e-14 trial_rc_gap=1e-14 old_div_gap=1e-15 trial_div_gap=1e-15 old_authoritative_phi_gap=1e-3 trial_authoritative_phi_gap=2e-3 same_state_operator=1 qualified=0 committed=0\n"
        c = "GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=1 merit=2 ratio=0.666 baseline=1 physical_converged=0 qualified=0\n"
        p = "GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=2 p_continuity=1e-2 k=1e-3 omega=1e-2 merit=2 converged=0 qualified=0\n"
        return t + c + p + "Converged NO Iterations 40\n"

    def test_valid_diagnostic(self):
        d = audit(self.fixture())
        self.assertEqual(d["iterations_committed"], 1)
        self.assertFalse(d["physical_qualification"])
        self.assertFalse(d["certified_convergence"])
        self.assertFalse(d["complete_40_iteration_run"])
        self.assertEqual(d["halted_trial"]["iteration"], 2)

    def test_forged_40_step_completion_rejected(self):
        self.assertFalse(audit(self.fixture())["complete_40_iteration_run"])

    def test_failed_iteration_proven(self):
        more = "GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=2 retry=0 beta=1 baseline=0 old_merit=2 trial_merit=1.99 ratio=0.9955 largest_growth=1 worst_block=0 joint_ok=1 monotone_ok=0 accepted=0 old_rc_gap=1e-14 trial_rc_gap=1e-14 old_div_gap=1e-15 trial_div_gap=1e-15 old_authoritative_phi_gap=1e-4 trial_authoritative_phi_gap=2e-4 same_state_operator=1 qualified=0 committed=0\n"
        d = audit(self.fixture() + more)
        self.assertEqual(d["halted_trial"]["iteration"], 2)
        self.assertTrue(d["halted_trial"]["all_rejected"])
        self.assertEqual(d["halted_trial"]["last_monotone_ok"], 0)

    def test_missing_trial(self):
        with self.assertRaises(ValueError):
            audit(self.fixture().split("\n", 1)[1])

    def test_fake_converged(self):
        with self.assertRaises(ValueError):
            audit(
                self.fixture().replace("physical_converged=0", "physical_converged=1")
            )

    def test_unclosed_rc(self):
        with self.assertRaises(ValueError):
            audit(self.fixture().replace("trial_rc_gap=1e-14", "trial_rc_gap=1e-5"))

    def test_claimed_qualification(self):
        with self.assertRaises(ValueError):
            audit(self.fixture().replace("qualified=0", "qualified=1"))

    def test_incompatible_trial_merit(self):
        with self.assertRaises(ValueError):
            audit(
                self.fixture().replace(
                    "merit=2 ratio=0.666 baseline", "merit=3 ratio=0.666 baseline"
                )
            )

    def test_unaccepted_commit(self):
        with self.assertRaises(ValueError):
            audit(self.fixture().replace("accepted=1", "accepted=0"))

    def test_broken_zero_step(self):
        with self.assertRaises(ValueError):
            audit(self.fixture().replace("old_div_gap=1e-15", "old_div_gap=1e-4"))


if __name__ == "__main__":
    unittest.main()
