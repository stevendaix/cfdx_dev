"""Mutation tests for fail-closed V57 6N Krylov + Schur audit."""
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from audit_v57_six_schur import audit_six, analyze_log, METHODS

MATRIX = 'COUPLED_MATRIX_FINGERPRINT_V40 A=abc b=def x0=ghi rows=608 nnz=7080 fnv64_not_sha256=1\n'
TRUTH = ('FPHYS_KRYLOV_TRUTH_V57 iteration=1 step=0 direction={d} '
         'arnoldi=0.93 true=0.93001 nonadditivity=0.001 gap=0.00001 '
         'consistency_limit=0.000001 qualified=0 production_committed=0\n')
SUMMARY = ('FPHYS_JFNK_V52 selected_schur={method} frozen_prestep=1 '
           'v53_global_block=1 central_jv=1 private_accepted=0 production_committed=0 '
           'qualified=0 v57_hybrid_jv_calls=33 v53_transport_nnz=1720 '
           'v57_true_krylov_checks=16 linear_true_relative=0.93001 '
           'v57_hybrid_consistent=0 v57_krylov_consistent=0 '
           'v57_hybrid_convection_gap=0.8 v57_hybrid_nonsmooth_faces=834 '
           'common_state_fingerprint=123 merit_before=218.05\n')
PCD=('COUPLED_LINEAR_UNCHANGED_RETRY_V37 iteration=1 schur=pcd committed=0 qualified=0\n'
     'CFDX production solver error: coupled momentum-continuity solve did not converge\n')

class TestV57Audit(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)
        for m in METHODS:
            contents=MATRIX
            if m=='pcd':contents+=PCD
            else:
                contents+=''.join(TRUTH.format(d=i) for i in range(1,17))
                contents+=SUMMARY.format(method=m)
                contents+='Converged NO Iterations 1\n'
            (self.folder/f'six_schur_{m}.log').write_text(contents)
    def alter(self,m,old,new):
        path=self.folder/f'six_schur_{m}.log'
        text=path.read_text()
        self.assertIn(old,text)
        path.write_text(text.replace(old,new,1))
    def must_fail(self):
        with self.assertRaises(ValueError):audit_six(self.folder)
    def test_clean_negative_campaign(self):
        r=audit_six(self.folder)
        self.assertFalse(r['qualified']);self.assertTrue(r['same_6n_initial_state'])
        self.assertEqual(len(r['methods']),6)
    def test_matrix_mismatch(self):
        self.alter('lsc','A=abc','A=other');self.must_fail()
    def test_state_mismatch(self):
        self.alter('bfbt','common_state_fingerprint=123','common_state_fingerprint=321');self.must_fail()
    def test_missing_direction(self):
        self.alter('simple',TRUTH.format(d=7),'');self.must_fail()
    def test_duplicate_direction(self):
        self.alter('simplec','direction=8 ','direction=7 ');self.must_fail()
    def test_wrong_krylov_budget(self):
        self.alter('lsc','v57_true_krylov_checks=16','v57_true_krylov_checks=15');self.must_fail()
    def test_false_qualification(self):
        self.alter('bfbt','qualified=0 production_committed=0','qualified=1 production_committed=0');self.must_fail()
    def test_false_newton_commit(self):
        self.alter('simple','private_accepted=0','private_accepted=1');self.must_fail()
    def test_bad_nonfinite_residual(self):
        self.alter('lsc','true=0.93001','true=nan');self.must_fail()
    def test_false_krylov_evidence(self):
        self.alter('simple','gap=0.00001','gap=0');self.must_fail()
    def test_bad_experimental_operator(self):
        self.alter('simplec','v57_hybrid_consistent=0','v57_hybrid_consistent=1');self.must_fail()
    def test_unqualified_physical_convergence(self):
        self.alter('lsc','Converged NO Iterations 1','Converged YES Iterations 1');self.must_fail()
    def test_pcd_hidden_fallback(self):
        self.alter('pcd','COUPLED_LINEAR_UNCHANGED_RETRY_V37','PCD_SILENT_FALLBACK');self.must_fail()
    def test_missing_fnv_label(self):
        self.alter('bfbt','fnv64_not_sha256=1','fnv64_not_sha256=0');self.must_fail()

if __name__=='__main__':unittest.main()
