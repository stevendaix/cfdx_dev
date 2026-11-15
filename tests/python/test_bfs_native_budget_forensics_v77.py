"""Fail-closed forensic contract for signed V75 budget / V76 private block."""
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from bfs_native_budget_forensics_v77 import check_native_log


def fixture(overrides=None):
    lines=[]
    for cycle in (1,2):
        for cell in (0,12):
            for eq in ('k','omega'):
                vals=(0.3,-0.2,-0.1,0.6,-0.5)
                raw=sum(vals)
                lines.append('BFS_SST_BUDGET_V75 iteration=23 retry=8 cycle='+str(cycle)+
                             ' cell='+str(cell)+' equation='+eq+
                             ' trace_selected=0 committed=0 '+
                             ' '.join(f'{nm}_delta={v:.17g}' for nm,v in zip(
                                 ('convection','diffusion','production','destruction','cross'),vals))+
                             f' terms_delta_sum={raw:.17g} raw_delta={raw:.17g} closure_error=0')
    for cycle,selected in ((1,1),(2,0)):
        lines.append(f'BFS_FIXED_DIRECTION_V74 iteration=23 retry=8 cycle={cycle} selected={selected}')
        lines.append('BFS_SST_BLOCK_V76 iteration=23 selected=0 reason=no_full_physical_descent')
    lines.append('BFS_FIXED_DIRECTION_SUMMARY_V74 iteration=23 native_momentum=0 native_sst=0 committed=0 qualified=0 mass_closed=0')
    for eq in ('k','omega'):
        lines.append(f'BFS_SST_ARGMAX_V75 iteration=23 cycle=2 equation={eq} after_cell=12')
    if overrides:
        for target,substitute in overrides:
            lines=[line.replace(target,substitute) for line in lines]
    return '\n'.join(lines)+'\n'


class V77NativeBudgetTests(unittest.TestCase):
    def run_fixture(self,content):
        with tempfile.TemporaryDirectory() as p:
            fname=Path(p)/'stderr.log'
            fname.write_text(content)
            return check_native_log(fname)

    def test_control(self):
        r=self.run_fixture(fixture())
        self.assertEqual(len(r['records']),8)
        self.assertFalse(r['native_momentum_certificate_for_candidate'])
        self.assertEqual(r['global_argmax_cycle2'],{'k':12,'omega':12})

    def test_reject_fabricated_certificates(self):
        with self.assertRaises(ValueError):
            self.run_fixture(fixture([('native_momentum=0','native_momentum=1')]))

    def test_reject_missing_budget(self):
        lines=fixture().splitlines()
        lines.pop(0)
        with self.assertRaises(ValueError):
            self.run_fixture('\n'.join(lines)+'\n')

    def test_reject_unbalanced_budget(self):
        with self.assertRaises(ValueError):
            self.run_fixture(fixture([('convection_delta=0.29999999999999999','convection_delta=0.79999999999999999')]))

    def test_reject_fabricated_sst_block_acceptance(self):
        with self.assertRaises(ValueError):
            self.run_fixture(fixture([('BFS_SST_BLOCK_V76 iteration=23 selected=0','BFS_SST_BLOCK_V76 iteration=23 selected=1')]))

    def test_reject_false_global_argmax(self):
        with self.assertRaises(ValueError):
            self.run_fixture(fixture([('after_cell=12','after_cell=0')]))


if __name__=='__main__':unittest.main()