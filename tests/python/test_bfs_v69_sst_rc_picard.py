import json
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from audit_bfs_v69_sst_rc_picard import audit

class V69AuditTest(unittest.TestCase):
    def test_guard_and_replay(self):
        with tempfile.TemporaryDirectory() as td:
            old = [f'GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=23 retry={i}' for i in range(149)]
            old += [f'BFS_SST_FLUX_SWAP_V68 iteration=23 retry={i}' for i in (0,8,10)]
            probe = ('BFS_SST_RC_PICARD_V69 iteration=23 retry=0 cycle=1 '
                     'beta=1 omega_phys_before=.5 omega_phys_after=.1 omega_ratio=.2 '
                     'k_phys_after=.01 merit_before=1 merit_after=.5 merit_ratio=.5 '
                     'continuity_before=.2 continuity_after=.2 rc_phi_change_linf=.03 '
                     'nut_change_linf=.01 rc_closed=1 v65_merit_ratio=.5 '
                     'v65_worst_growth=.6 v65_worst_block=0 v65_joint_ok=1 '
                     'v65_monotone_ok=1 v65_would_pass=1 '
                     'private_transport=1 committed=0 qualified=0')
            base = Path(td)/'base.log'; final = Path(td)/'new.log'
            base.write_text('\n'.join(old)+'\n'); final.write_text('\n'.join(old+[probe])+'\n')
            result = audit(final,base)
            self.assertEqual(len(result['probe_records']),1)
            self.assertFalse(result['qualified'])
            final.write_text('\n'.join(old+[probe.replace('committed=0','committed=1')])+'\n')
            with self.assertRaises(ValueError): audit(final,base)
            final.write_text('\n'.join(old[:-1]+[probe])+'\n')
            with self.assertRaises(ValueError): audit(final,base)
            final.write_text('\n'.join(old+[probe.replace('rc_closed=1','rc_closed=0')])+'\n')
            with self.assertRaises(ValueError): audit(final,base)
if __name__ == '__main__': unittest.main()