"""Audit actual V13 committed global rows against linear/SST statuses.

No renormalization or retuning. Detect non-monotone physical momentum,
report the best/final residual, and explicitly deny a convergence claim
unless the solver itself supplies a qualified final state. Read-only.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path


def fields(line: str) -> dict[str,str]:
    return dict(re.findall(r'([A-Za-z][A-Za-z_0-9]*)=([^\s]+)',line))

def audit(text: str) -> dict:
    records={}
    for line in text.splitlines():
        if not line.startswith('GLOBAL_NONLINEAR_COMMITTED_V13 '):
            continue
        d=fields(line)
        required={'iteration','U','p_continuity','k','omega','merit','converged','qualified'}
        if required-d.keys():raise ValueError('truncated V13 committed row')
        i=int(d['iteration']);
        if i<1 or i in records:raise ValueError('nonpositive or duplicate committed iteration')
        row={'iteration':i}
        for x in ('U','p_continuity','k','omega','merit'):
            v=float(d[x]);
            if not math.isfinite(v) or v<0:raise ValueError('invalid committed physical '+x)
            row[x]=v
        for x in ('converged','qualified'):
            if d[x] not in ('0','1'):raise ValueError('invalid V13 certification flag')
            row[x]=d[x]=='1'
        records[i]=row
    if not records:raise ValueError('no committed V13 physical evidence')
    if sorted(records)!=list(range(1,len(records)+1)):
        raise ValueError('missing or out-of-order committed V13 iteration')
    rows=[records[i] for i in sorted(records)]
    worsening=[row['iteration'] for last,row in zip(rows,rows[1:])
               if row['U']>last['U']]
    if any(row['qualified'] and not row['converged'] for row in rows):
        raise ValueError('inconsistent V13 qualification flags')
    lin=[];sst=[]; post=[]
    for line in text.splitlines():
        if line.startswith('BFS_PHASE_V10 ') and 'phase=COUPLED_LINEAR_END' in line:
            d=fields(line)
            if {'iteration','true_relative','status'}<=d.keys():
                lin.append((int(d['iteration']),float(d['true_relative']),d['status']))
        elif line.startswith('BFS_PHASE_V10 ') and 'phase=SST_TRANSPORT_END' in line:
            d=fields(line)
            if {'iteration','converged','k_residual','omega_residual'}<=d.keys():
                sst.append((int(d['iteration']),d['converged'],float(d['k_residual']),
                            float(d['omega_residual'])))
        elif line.startswith('POST_SST_MOMENTUM '):
            d=fields(line)
            if {'iteration','relative_nut_change','old_nut_operator',
                'refreshed_nut_operator'}<=d.keys():
                post.append((int(d['iteration']),float(d['relative_nut_change']),
                             float(d['old_nut_operator']),float(d['refreshed_nut_operator'])))
    if len(lin)!=len(rows) or len(sst)!=len(rows) or len(post)!=len(rows):
        raise ValueError('missing matched linear, SST or post-SST reassembly records')
    for group in (lin,sst,post):
        if [r[0] for r in group]!=list(range(1,len(rows)+1)):
            raise ValueError('iteration mismatch across coupled/SST/V13 audit')
        if any(any(not math.isfinite(float(x)) for x in r[1:] if x not in ('0','1')) for r in group):
            raise ValueError('nonfinite supporting linear/SST record')
    best=min(rows,key=lambda r:r['U'])
    return {
        'schema':'CFDX-BFS-COMMITTED-CONTRACTION/V64',
        'status':'OBSERVED_PRODUCTION_ITERATIONS_NOT_INDEPENDENT_PHYSICAL_VALIDATION',
        'count':len(rows),
        'momentum_increase_iterations':worsening,
        'nonmonotone_increases':len(worsening),
        'starting_momentum':rows[0]['U'],
        'final_momentum':rows[-1]['U'],
        'best_momentum':best['U'],
        'best_iteration':best['iteration'],
        'momentum_drop_factor':rows[0]['U']/rows[-1]['U'] if rows[-1]['U'] else None,
        'final_continuity':rows[-1]['p_continuity'],
        'max_linear_true_relative':max(x[1] for x in lin),
        'max_sst_k_residual':max(x[2] for x in sst),
        'max_sst_omega_residual':max(x[3] for x in sst),
        'all_linear_steps_status_zero':all(x[2]=='0' for x in lin),
        'all_sst_steps_converged':all(x[1]=='1' for x in sst),
        'last_relative_nut_change':post[-1][1],
        'solver_qualified_last':rows[-1]['qualified'],
        'solver_converged_last':rows[-1]['converged'],
        'iterations':rows,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    opts=p.parse_args()
    report=audit(opts.log.read_text(errors='replace'))
    opts.output.parent.mkdir(parents=True,exist_ok=True)
    opts.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(f"V64_COMMITTED_DIAG count={report['count']} "
          f"rises={report['nonmonotone_increases']} "
          f"last_qualified={int(report['solver_qualified_last'])}")

if __name__=='__main__':main()