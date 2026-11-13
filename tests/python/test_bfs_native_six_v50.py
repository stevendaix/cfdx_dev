"""Mutation tests of the V50 fail-closed forensic contract."""
import pathlib
import sys
import unittest
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from analyze_bfs_native_six_v50 import verify
T='NATIVE_SIX_VECTOR_V50 iteration=1 retry=0 closed=1 cells=2 entries=12 U=1.25 continuity=0.008 k=0.06 omega=0.05 merit=1.25 worst_u=0 worst_v=1 worst_k=0 worst_omega=1 fingerprint=1234 parity=1 same_native_flux=1 experimental=1 committed=0 qualified=0'
N='NATIVE_RC_ORACLE_V49 iteration=1 retry=0 closed=1 fixed_point_iterations=25 fixed_point_defect=1e-12 repeated_same_state=1 native_flux_scale_linf=0.001 U=1.25 continuity=0.008 k=0.06 omega=0.05 carrier_flux_gap_l1=0.03 carrier_flux_gap_linf=0.001 worst_face=1 full_six_field_oracle_qualified=0 physical_converged=0 committed=0 qualified=0'
C='NATIVE_SIX_COMMITTED_V50 iteration=1 cells=2 entries=12 fingerprint=1234 merit=1.25 U=1.25 continuity=0.008 k=0.06 omega=0.05 flux_gap_l1=0.03 committed=1 physical_converged=0 qualified=0'
G='NATIVE_SIX_FINAL_GATE_V50 iteration=1 native_available=1 eligible=0 gap=0.03 tolerance=1e-8 experimental=1 qualified=0'
V='GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=1.9 p_continuity=0.01 k=1e-12 omega=1e-12 nut_closure=0 merit=1.9 converged=0 qualified=0'
LOG='\n'.join((T,N,C,V))
class Forensics(unittest.TestCase):
    def test_valid(self):
        r=verify(LOG)
        self.assertEqual(r['status'],'DIAGNOSTIC_ONLY_PHYSICAL_UNQUALIFIED')
        self.assertFalse(r['qualified'])
        self.assertEqual(r['committed_outer_states'],1)
    def rejects(self, bad):
        with self.assertRaises(ValueError):verify(bad)
    def test_final_gate_valid(self):
        self.assertEqual(verify(LOG+'\n'+G)['committed_outer_states'],1)
    def test_final_gate_false_accept(self):
        self.rejects(LOG+'\n'+G.replace('eligible=0','eligible=1'))
    def test_final_gate_wrong_gap(self):
        self.rejects(LOG+'\n'+G.replace('gap=0.03','gap=0.04'))
    def test_missing_commit(self):self.rejects(LOG.replace(C,''))
    def test_false_qualification(self):self.rejects(LOG.replace('qualified=0','qualified=1'))
    def test_fake_production_converged(self):
        self.rejects(LOG.replace('converged=0 qualified=0','converged=1 qualified=0'))
    def test_mismatched_fingerprint(self):self.rejects(LOG.replace('fingerprint=1234 merit','fingerprint=1235 merit'))
    def test_bad_shape(self):self.rejects(LOG.replace('entries=12','entries=10'))
    def test_missing_production(self):self.rejects(LOG.replace(V,''))
    def test_absent_v49(self):self.rejects(LOG.replace(N,''))
    def test_bad_merit(self):self.rejects(LOG.replace('merit=1.25 worst','merit=0.05 worst'))
    def test_k_mismatch(self):self.rejects(LOG.replace('k=0.06 omega=0.05 merit','k=0.16 omega=0.05 merit'))
    def test_nonfinite(self):self.rejects(LOG.replace('continuity=0.008','continuity=nan'))
    def test_duplicate(self):self.rejects(LOG+'\n'+T)
    def test_closed_false(self):self.rejects(LOG.replace('closed=1','closed=0'))
    def test_no_parity(self):self.rejects(LOG.replace('parity=1','parity=0'))
    def test_same_flux_broken(self):self.rejects(LOG.replace('same_native_flux=1','same_native_flux=0'))
    def test_v49_false_physical(self):self.rejects(LOG.replace('physical_converged=0','physical_converged=1'))
if __name__=='__main__':unittest.main()
