import unittest
from scripts.audit_bfs_causal_v66 import audit, check_unmodified_v65_trace

class BfsCausalV66Tests(unittest.TestCase):
    def sample(self):
        v65='GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=23 retry=0 beta=1 old_merit=2 trial_merit=1 ratio=0.5 largest_growth=20 worst_block=5 joint_ok=0 monotone_ok=1 accepted=0 old_authoritative_phi_gap=0.002 trial_authoritative_phi_gap=0.003 committed=0 qualified=0\n'
        car='BFS_CARRIER_COHERENCE_V66 iteration=23 retry=0 beta=1 old_carrier_gap=0.002 trial_carrier_gap=0.003 old_carrier_div_gap=0.0001 trial_carrier_div_gap=0.0002 worst_face=5 worst_div_cell=7 old_rau_carrier_gap=0.01 trial_rau_carrier_gap=0.3 old_rau_carrier_available=1 old_nut_carrier_gap=0.02 trial_nut_carrier_gap=0.04 rau_change=0.2 rau_cell=3 rau_component=1 nut_change=0.4 nut_cell=8 omega_change_cell=7 omega_trial_cell=7 omega_old_cell=6 physical_qualified=0 committed=0\n'
        def hot(c):
            return f'BFS_OMEGA_HOTSPOT_V66 iteration=23 retry=0 beta=1 cell={c} x=0.2 y=0.3 z=0 old_raw=1 trial_raw=2 delta=1 omega_old=1 omega_trial=2 k_old=0.3 k_trial=0.4 nut_old=0.1 nut_trial=0.2 physical_qualified=0 committed=0\n'
        bud='BFS_OMEGA_BUDGET_V66 iteration=23 retry=0 beta=1 cell=7 convection=1 diffusion=1 production=-2 destruction=1 cross=0 sum=1 residual_delta=1 closure=0 old_row=1 trial_row=2 old_wall_flux=0 trial_wall_flux=0 owner_constraint=0 oracle_gap_old=0 oracle_gap_trial=0 physical_qualified=0 committed=0\n'
        return v65+car+hot(7)+hot(7)+hot(6)+bud
    def test_valid(self):
        r=audit(self.sample());self.assertEqual(r['trial_count'],1)
        self.assertEqual(r['budget_count'],1);self.assertFalse(r['physical_qualification'])
    def test_negative_roundoff_closure(self):
        result=audit(self.sample().replace("closure=0", "closure=-1e-14"))
        self.assertEqual(result["budget_count"], 1)
    def test_significant_negative_closure_rejected(self):
        self.assertRaises(ValueError, audit, self.sample().replace("closure=0", "closure=-1e-2"))
    def test_missing(self):
        t=self.sample();self.assertRaises(ValueError,audit,t.replace('BFS_OMEGA_BUDGET_V66','BAD_OMEGA_BUDGET_V66'))
    def test_inconsistent_flux(self):
        self.assertRaises(ValueError,audit,self.sample().replace('trial_carrier_gap=0.003','trial_carrier_gap=0.004'))
    def test_fake_omega(self):
        self.assertRaises(ValueError,audit,self.sample().replace('delta=1','delta=3'))
    def test_wrong_residual(self):
        self.assertRaises(ValueError,audit,self.sample().replace('residual_delta=1','residual_delta=2'))
    def test_forged_pass(self):
        self.assertRaises(ValueError,audit,self.sample().replace('committed=0 qualified=0','committed=1 qualified=1'))
    def test_missing_cell(self):
        self.assertRaises(ValueError,audit,self.sample().replace('cell=7','cell=153'))
    def test_duplicate(self):
        self.assertRaises(ValueError,audit,self.sample()+self.sample())
    def test_nonfinite(self):
        self.assertRaises(ValueError,audit,self.sample().replace('rau_change=0.2','rau_change=nan'))
    def test_reference_parity(self):
        text=self.sample()+'GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=1 merit=1\n'
        one=text+'GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 merit=1\n'
        two=text+'GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 merit=2\n'
        self.assertTrue(check_unmodified_v65_trace(one,one))
        self.assertRaises(ValueError,check_unmodified_v65_trace,one,two)
    def test_missing_iteration23(self):
        self.assertRaises(ValueError,audit,self.sample().replace('iteration=23','iteration=22'))

if __name__=='__main__': unittest.main()