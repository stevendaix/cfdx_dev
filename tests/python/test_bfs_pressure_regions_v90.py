@@ -0,0 +1,20 @@
import tempfile, unittest
from pathlib import Path
from scripts.audit_bfs_pressure_regions_v90 import accepted_iterations
class TestV90PressureAudit(unittest.TestCase):
    def write(self,text):
        tmp=tempfile.TemporaryDirectory();p=Path(tmp.name)/'log';p.write_text(text)
        try:return accepted_iterations(p)
        finally:tmp.cleanup()
    def test_accept_requires_native_commits(self):
        self.assertEqual(self.write(
            'GLOBAL_PROGRESS_TRIAL_V85 iteration=1 accepted=1\n'
            'GLOBAL_PROGRESS_COMMITTED_V85 iteration=1 committed=1 converged=0 qualified=0\n'
            'GLOBAL_PROGRESS_COMMITTED_V85 iteration=2 committed=2 converged=0 qualified=0\n'),[1,2])
    def test_no_trial_counted(self):
        with self.assertRaises(ValueError):self.write('GLOBAL_PROGRESS_TRIAL_V85 iteration=1 accepted=1\n')
    def test_reject_nonmonotone(self):
        with self.assertRaises(ValueError):self.write('GLOBAL_PROGRESS_COMMITTED_V85 committed=1 converged=0 qualified=0\nGLOBAL_PROGRESS_COMMITTED_V85 committed=1 converged=0 qualified=0\n')
    def test_reject_certification(self):
        with self.assertRaises(ValueError):self.write('GLOBAL_PROGRESS_COMMITTED_V85 committed=1 converged=1 qualified=1\n')
if __name__=='__main__':unittest.main()
