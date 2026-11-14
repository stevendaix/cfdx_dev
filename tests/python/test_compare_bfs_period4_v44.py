import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
from compare_bfs_period4_v44 import compare

FP="COUPLED_MATRIX_FINGERPRINT_V40 A=aaa b=bbb x0=ccc rows=608 nnz=7000\n"
def commit(n,merit):
    return (f"GLOBAL_NONLINEAR_COMMITTED_V13 iteration={n} U={merit} "
            f"p_continuity=0.001 k=1e-11 omega=1e-13 nut_closure=0 "
            f"merit={merit} converged=0 qualified=0\n")
def trace(n,phase):
    return f"OUTER_FIELD_PERIOD_V43 iteration={n} period_four_candidate={phase}\n"
TRIGGER=("COUPLED_PERIOD4_FULL_V44 iteration=6 retry=0 k_omega_averaged=1 "
         "nut_recomputed=1 strain_recomputed=1 pressure_gauge_preserved=1 "
         "committed=0 physical_qualified=0\n")
def log(changed=False,with_trigger=True,corrupt=None):
    text=FP
    for i in range(1,9):
        if with_trigger and changed and i==6: text+=TRIGGER
        text+=commit(i, (9.0-i)*1e-2 + (0.002 if changed and i>=6 else 0.0))
        if i>=5: text+=trace(i,1 if i==5 else 0)
    text+="Converged NO Iterations 8\n"
    return text if corrupt is None else text.replace(*corrupt)

class TestV44Comparison(unittest.TestCase):
    def test_matched_preintervention_and_true_return(self):
        data=compare(log(),log(True),True)
        self.assertEqual(data['first_full_sst_predictor'],6)
        self.assertEqual(data['common_committed_states'],8)
        self.assertFalse(data['physical_converged'])
        self.assertTrue(data['merit_ratio_at_first_trial']>1)
    def test_fail_closed_wrong_matrix(self):
        with self.assertRaisesRegex(ValueError,"fingerprint"):
            compare(log(), log(True).replace("A=aaa","A=xxx"),True)
    def test_fail_closed_without_trigger(self):
        with self.assertRaisesRegex(ValueError,"never triggered"):
            compare(log(),log(),True)
    def test_fail_closed_fabricated_convergence(self):
        with self.assertRaisesRegex(ValueError,"certification"):
            compare(log(),log(True).replace("qualified=0","qualified=1"),True)
    def test_fail_closed_noncommitted_predictor(self):
        with self.assertRaisesRegex(ValueError,"precedes|preceding"):
            compare(log(),log(True).replace("period_four_candidate=1","period_four_candidate=0"),True)
    def test_fail_closed_drift_before_trigger(self):
        with self.assertRaisesRegex(ValueError,"BEFORE"):
            compare(log(),log(True).replace("iteration=4 U=0.05","iteration=4 U=0.5").replace("merit=0.05","merit=0.5"),True)

if __name__=='__main__': unittest.main()
