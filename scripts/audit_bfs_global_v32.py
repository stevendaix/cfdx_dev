#!/usr/bin/env python3
"""Fail-closed external nonlinear convergence audit for CFDX BFS V32.

An SST internal accepted step never counts as an external U,p,k,omega commit.
The only evidence counted is GLOBAL_NONLINEAR_COMMITTED_V13, and every
committed vector must be bounded/finite and monotonically contract on the
same final operator. This is NOT a substitute for a NASA validation oracle.
"""
from __future__ import annotations
import argparse
import json
import math
import re
import sys
from pathlib import Path

_PATTERN=re.compile(r'^GLOBAL_NONLINEAR_COMMITTED_V13\s+(.+)$',re.M)
_PAIRS=re.compile(r'([A-Za-z_]+)=([^ ]+)')
_REQUIRED=('iteration','U','p_continuity','k','omega','nut_closure','merit','converged')

def audit(log_text:str, min_commits:int=3, tolerance:float=1.e-8):
    if not isinstance(min_commits,int) or min_commits<1:
        raise ValueError('invalid required number of commits')
    if not math.isfinite(tolerance) or tolerance<=0:
        raise ValueError('invalid physical threshold')
    commits=[]
    for match in _PATTERN.finditer(log_text):
        row=dict(_PAIRS.findall(match.group(1)))
        if any(k not in row for k in _REQUIRED):
            raise ValueError('incomplete committed-state evidence')
        try:
            fields={k:float(row[k]) for k in _REQUIRED if k not in ('iteration','converged')}
            iteration=int(row['iteration'])
            converged=int(row['converged'])
        except (ValueError,OverflowError) as exc:
            raise ValueError('invalid committed-state numeric data') from exc
        if iteration<=0 or converged not in (0,1) or any(
            not math.isfinite(x) or x<0 for x in fields.values()):
            raise ValueError('invalid committed-state physical metrics')
        physical_max=max(fields[k] for k in ('U','p_continuity','k','omega','nut_closure'))
        if not math.isclose(physical_max,fields['merit'],rel_tol=1.e-6,abs_tol=1.e-12):
            raise ValueError('claimed global merit differs from physical components')
        if converged and physical_max>tolerance:
            raise ValueError('false physical convergence certificate')
        if commits:
            if iteration<=commits[-1]['iteration']:
                raise ValueError('non-increasing committed outer iteration')
            if not (fields['merit']<commits[-1]['merit']) and not converged:
                raise ValueError('committed physical residual did not contract')
        commits.append({'iteration':iteration,**fields,'converged':bool(converged)})
    # A successful software run with accepted outer states is still not a
    # *physically converged* steady BFS until all independent gates certify.
    reached_target=(len(commits)>=min_commits and commits[-1]['converged']
                    and all(commits[-1][k]<=tolerance for k in
                            ('U','p_continuity','k','omega','nut_closure')))
    return {'status':'PHYSICALLY_CONVERGED' if reached_target else 'INCOMPLETE_UNCONVERGED',
            'qualified':reached_target,'committed_outer_states':len(commits),
            'required_commits':min_commits,'last_committed':commits[-1] if commits else None,
            'sst_internal_accepts':len(re.findall(r'^SST_JFNK_V22[^\n]*accepted=1',log_text,re.M)),
            'outer_rejections':len(re.findall(r'^CFDX nonlinear iteration \d+ rejected;',log_text,re.M))}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('logfile',type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--minimum-commits',type=int,default=3)
    p.add_argument('--tolerance',type=float,default=1.e-8)
    a=p.parse_args()
    try:
        result=audit(a.logfile.read_text(encoding='utf-8',errors='replace'),a.minimum_commits,a.tolerance)
    except (OSError,ValueError) as e:
        result={'status':'FAIL_INVALID_EVIDENCE','qualified':False,'reason':str(e)}
    data=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if a.output:a.output.write_text(data,encoding='utf-8')
    print(data,end='')
    return 0 if result['qualified'] else 3
if __name__=='__main__':sys.exit(main())
