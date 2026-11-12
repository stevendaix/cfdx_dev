import importlib.util
from pathlib import Path
import sys
import unittest

SCRIPT=Path(__file__).resolve().parents[2]/"scripts"/"analyze_bfs_face_forensics_v46.py"
spec=importlib.util.spec_from_file_location("bfs_v46_face",SCRIPT)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)

ROW=("POST_SST_FACE_CELL_V46 iteration=59 component=1 cell=2 faces=1 "
     "u=2 bounded_mode=1 conv=4 visc=0 transpose=0 pressure=1 body=-0.25 "
     "mass_div=2 bounded=-4 collocated_delta=2 physical=4.75 "
     "face_operator=0.75 csr=0.75 closure=0 closure_scale=10 "
     "face_ok=1 csr_ok=1 all_csr_ok=1 nu_t=0.005 k=3 omega=500 "
     "committed=0 qualified=0")
FACE=("POST_SST_FACE_V46 iteration=59 component=1 cell=2 face=7 owner=2 "
      "neighbour=-1 sign=1 patch=wall mass=2 conv=4 visc=0 transpose=0 "
      "collocated_delta=2 committed=0 qualified=0")

class EvidenceTests(unittest.TestCase):
    def test_valid_evidence_only(self):
        d=mod.certify(ROW+'\n'+FACE)
        self.assertEqual(d['status'],'DIAGNOSTIC_ONLY_UNQUALIFIED')
        self.assertFalse(d['physically_converged'])
        self.assertEqual(d['faces'],1)
        self.assertEqual(d['iterations'],[59])
    def test_consecutive_phase_decomposition(self):
        next_row=(ROW.replace('iteration=59','iteration=60')
                  .replace('conv=4 ', 'conv=5 ')
                  .replace('physical=4.75 ', 'physical=5.75 ')
                  .replace('face_operator=0.75 ', 'face_operator=1.75 ')
                  .replace('csr=0.75 ', 'csr=1.75 '))
        next_face=FACE.replace('iteration=59','iteration=60').replace('conv=4 ', 'conv=5 ')
        d=mod.certify(ROW+'\n'+FACE+'\n'+next_row+'\n'+next_face)
        change=d['consecutive_phase_deltas']
        self.assertEqual(len(change),1)
        self.assertEqual((change[0]['from'],change[0]['to']),(59,60))
        self.assertAlmostEqual(change[0]['delta_csr'],1.0)
        self.assertAlmostEqual(change[0]['delta_terms']['conv'],1.0)
        self.assertEqual(change[0]['dominant_terms'][0],'conv')
        self.assertAlmostEqual(change[0]['most_changed_faces'][0]['delta_fv_face'],1.0)
    def test_mass_flux_direction_switch(self):
        # Preserve the bounded term and the convective ledger while flipping
        # the signed oriented mass flux (zero velocity coefficient here).
        row_a=ROW.replace('u=2 bounded_mode=1','u=0 bounded_mode=1')
        row_a=row_a.replace('bounded=-4 ','bounded=0 ').replace('face_operator=0.75 ','face_operator=4.75 ').replace('csr=0.75 ','csr=4.75 ')
        face_a=FACE
        row_b=row_a.replace('iteration=59','iteration=60').replace('mass_div=2 ','mass_div=-2 ')
        face_b=FACE.replace('iteration=59','iteration=60').replace('mass=2 ','mass=-2 ')
        result=mod.certify(row_a+'\n'+face_a+'\n'+row_b+'\n'+face_b)
        self.assertEqual(result['consecutive_phase_deltas'][0]['mass_flux_direction_reversals'],1)
    def test_phase_gap_does_not_pretend_continuity(self):
        later_row=ROW.replace('iteration=59','iteration=61')
        later_face=FACE.replace('iteration=59','iteration=61')
        d=mod.certify(ROW+'\n'+FACE+'\n'+later_row+'\n'+later_face)
        self.assertFalse(d['consecutive_phase_deltas'])
    def test_changed_face_identity_rejected(self):
        later_row=ROW.replace('iteration=59','iteration=60')
        later_face=FACE.replace('iteration=59','iteration=60').replace('patch=wall','patch=step_wall')
        with self.assertRaises(ValueError):
            mod.certify(ROW+'\n'+FACE+'\n'+later_row+'\n'+later_face)
    def test_nonfinite(self):
        with self.assertRaises(ValueError):
            mod.certify(ROW.replace('nu_t=0.005','nu_t=nan')+'\n'+FACE)
    def test_false_certificate(self):
        with self.assertRaises(ValueError):
            mod.certify(ROW.replace('qualified=0','qualified=1')+'\n'+FACE)
    def test_missing_face(self):
        with self.assertRaises(ValueError): mod.certify(ROW)
    def test_duplicate_face(self):
        with self.assertRaises(ValueError): mod.certify(ROW+'\n'+FACE+'\n'+FACE)
    def test_bad_mass(self):
        with self.assertRaises(ValueError):
            mod.certify(ROW+'\n'+FACE.replace('mass=2','mass=3'))
    def test_bad_pressure(self):
        with self.assertRaises(ValueError):
            mod.certify(ROW.replace('pressure=1','pressure=2')+'\n'+FACE)
    def test_bad_orientation(self):
        with self.assertRaises(ValueError):
            mod.certify(ROW+'\n'+FACE.replace('sign=1','sign=-1'))
    def test_bad_closure(self):
        with self.assertRaises(ValueError):
            mod.certify(ROW.replace('csr=0.75','csr=0.8')+'\n'+FACE)
    def test_explicit_diagnostic_error(self):
        with self.assertRaises(ValueError):
            mod.certify(ROW+'\n'+FACE+'\nPOST_SST_FACE_ERROR_V46 iteration=59 reason=bad qualified=0')

if __name__ == '__main__': unittest.main()
