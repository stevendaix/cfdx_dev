from __future__ import annotations
import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from qualify_bfs_two_grids_v40 import compile_report, read_pcd_csv
from generate_th03_bfs_case import SCHUR_MODELS


def fake(cells):
    cases={x:{'n_cells':cells,'n_faces':cells*3,'physical_sha256':'same'} for x in SCHUR_MODELS}
    results=[{'schur':x,'execution':'TIMEOUT_INCOMPLETE','nonlinear_convergence_proven':False,
              'operator_fingerprint_v40':None,'first_blocker':{'first_blocker':'KRYLOV_SCHUR_RUNNING_TIMEOUT'},
              'wall_seconds':1.0} for x in SCHUR_MODELS]
    return {'identical_physical_inputs':True,'cases':cases,'results':results,
            'all_six_identical_original_operator':False,'source_kind':'GENERATED'}

class TestBFSV40(unittest.TestCase):
    def test_timeout_never_passes(self):
        result=compile_report(fake(712),fake(20540))
        self.assertFalse(result['all_twelve_outer_converged'])
        self.assertFalse(result['nasa_validation_qualified'])
    def test_malformed_success_fails_closed(self):
        a=fake(712);a['results'][0]['nonlinear_convergence_proven']=True
        with self.assertRaisesRegex(ValueError,'invalid converged claim'):
            compile_report(a,fake(20540))
    def test_missing_schur_fails_closed(self):
        a=fake(712);a['results'].pop()
        with self.assertRaisesRegex(ValueError,'missing or duplicate'):
            compile_report(a,fake(20540))
    def test_unmatched_case_fails_closed(self):
        a=fake(712);a['cases']['pcd']['physical_sha256']='changed'
        with self.assertRaisesRegex(ValueError,'physical input mismatch'):
            compile_report(a,fake(20540))
if __name__=='__main__':unittest.main()
