from __future__ import annotations
import sys
from pathlib import Path
import unittest
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from bfs_schur_fingerprint_v40 import parse_fingerprint,compare_six_fingerprints
from generate_th03_bfs_case import SCHUR_MODELS

LINE = ('COUPLED_MATRIX_FINGERPRINT_V40 A=0123456789abcdef b=0123456789abcdef '
        'x0=0123456789abcdef rows=608 nnz=12345 fnv64_not_sha256=1 '
        'physical_converged=0')

class TestFingerprint(unittest.TestCase):
    def test_capture_parse(self):
        self.assertEqual(parse_fingerprint(LINE)['rows'],608)
        self.assertIsNone(parse_fingerprint('not a capture'))
    def test_malformed_rejected(self):
        self.assertIsNone(parse_fingerprint(LINE.replace('fnv64_not_sha256=1','fnv64_not_sha256=0')))
        with self.assertRaises(ValueError): parse_fingerprint(LINE+'\n'+LINE)
    def test_common_matrix(self):
        rows=[{'schur':key,'operator_fingerprint_v40':parse_fingerprint(LINE)} for key in SCHUR_MODELS]
        result=compare_six_fingerprints(rows,SCHUR_MODELS)
        self.assertTrue(result['all_six_same_A_b_x0_fingerprint'])
        self.assertFalse(result['identity_is_cryptographic_proof'])
    def test_missing_or_changed_fails_closed(self):
        rows=[{'schur':key,'operator_fingerprint_v40':parse_fingerprint(LINE)} for key in SCHUR_MODELS]
        rows[-1]['operator_fingerprint_v40']=None
        self.assertFalse(compare_six_fingerprints(rows,SCHUR_MODELS)['all_six_same_A_b_x0_fingerprint'])
        rows[-1]['operator_fingerprint_v40']=parse_fingerprint(LINE.replace('rows=608','rows=609'))
        self.assertFalse(compare_six_fingerprints(rows,SCHUR_MODELS)['all_six_same_A_b_x0_fingerprint'])
if __name__=='__main__': unittest.main()
