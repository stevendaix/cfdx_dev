@@ -0,0 +1,39 @@
import tempfile
import unittest
from pathlib import Path
from scripts.plot_bfs_sst_v82 import parse_private_cycle


class TestEvidenceV82(unittest.TestCase):
    def make(self, *, v82=True, alter=None):
        body='''PRESSURE_HILO_PHYSICAL_FLUX stored_fp64_relative=8.87e-14
INCOMPRESSIBLE_ITER momentum_residual_rhs_relative=0.58
BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_BEGIN budget=10 target=1e-8
SST_TRANSPORT_DIAGNOSTIC iteration=1 k_returned_row_scaled=0.05 omega_returned_row_scaled=800 k_local_relative=0.4 omega_local_relative=1.2
SST_IMPLICIT_V82 iteration=1 selected=picard_full accepted=1 alpha=1 evaluations=4 k_before=0.7 k_after=0.8 omega_before=400 omega_after=1 old_merit=400 new_merit=1 reason=reassembled_global_max_descent fixed_rhie_chow=1 outer_U_p_flux_qualified=0 physical_qualified=0
BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_END committed=0 converged=0 iterations=10 k_residual=0.01 omega_residual=40
'''
        if not v82:body='\n'.join(l for l in body.splitlines() if not l.startswith('SST_IMPLICIT_V82'))+'\n'
        if alter:body=alter(body)
        f=Path(self.tmp.name)/('v82.log' if v82 else 'baseline.log');f.write_text(body);return f
    def setUp(self):self.tmp=tempfile.TemporaryDirectory()
    def tearDown(self):self.tmp.cleanup()
    def test_synthetic_valid_but_not_physical(self):
        r=parse_private_cycle(self.make(),implicit=True)
        self.assertEqual(r['sst_committed'],0)
        self.assertEqual(r['updates'][0]['selected'],'picard_full')
        self.assertEqual(parse_private_cycle(self.make(v82=False))['sst_iterations'],10)
    def reject(self,alter):
        with self.assertRaises(ValueError):parse_private_cycle(self.make(alter=alter),implicit=True)
    def test_no_false_commit(self):self.reject(lambda s:s.replace('committed=0','committed=1'))
    def test_no_lax_threshold(self):self.reject(lambda s:s.replace('target=1e-8','target=1e-5'))
    def test_no_pressure_failure(self):self.reject(lambda s:s.replace('8.87e-14','1e-4'))
    def test_no_missing_end(self):self.reject(lambda s:s.replace('phase=SST_TRANSPORT_END','phase=NEVER'))
    def test_no_false_candidate_certificate(self):self.reject(lambda s:s.replace('physical_qualified=0','physical_qualified=1'))
    def test_no_false_merit_descent(self):self.reject(lambda s:s.replace('new_merit=1 ','new_merit=500 '))
    def test_no_nonfinite(self):self.reject(lambda s:s.replace('omega_local_relative=1.2','omega_local_relative=nan'))
    def test_no_wrong_budget(self):self.reject(lambda s:s.replace('budget=10','budget=200'))
    def test_no_missing_v82_evidence(self):
        with self.assertRaises(ValueError):parse_private_cycle(self.make(v82=False),implicit=True)

if __name__=='__main__':unittest.main()
