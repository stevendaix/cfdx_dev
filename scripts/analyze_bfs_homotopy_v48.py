"""Fail-closed forensic verifier for the experimental V48 Rhie-Chow path.

No SST internal certificate or V13/V14 commit by itself qualifies Driver-Seegmiller.
This parser deliberately never promotes a path-dependent homotopy result to
physical qualification without a separate state-only F(X) oracle.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re

NUMBER = r"([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)"
TRIES = re.compile(r"^COUPLED_RAU_HOMOTOPY_V48 (?P<items>.+)$",re.M)
NATIVE = re.compile(r"^COUPLED_NATIVE_FLUX_GAP_V48 (?P<items>.+)$",re.M)
COMMITS = re.compile(r"^GLOBAL_NONLINEAR_COMMITTED_V13 (?P<items>.+)$",re.M)
TRIALS = re.compile(r"^COUPLED_COMPOSITE_TRIAL_V47 (?P<items>.+)$",re.M)

def fields(items: str) -> dict[str,str]:
    return dict(re.findall(r"([a-zA-Z][\w]*)=([^\s]+)",items))

def numeric(data:dict[str,str],key:str) -> float:
    if key not in data:
        raise ValueError(f"missing {key}")
    v=float(data[key])
    if not math.isfinite(v):
        raise ValueError(f"non-finite {key}")
    return v

def inspect(log: str, min_commits:int=1) -> dict:
    attempt_rows=[fields(m['items']) for m in TRIES.finditer(log)]
    trial_rows=[fields(m['items']) for m in TRIALS.finditer(log)]
    commit_rows=[fields(m['items']) for m in COMMITS.finditer(log)]
    native_rows=[fields(m['items']) for m in NATIVE.finditer(log)]
    if not attempt_rows or not commit_rows or not trial_rows:
        raise ValueError("missing V48 trials, V47 carrier trials or committed V13 evidence")
    # First outer state only establishes a committed rAU reference: there is
    # deliberately no prior coefficient to continue. From state 2 on, every
    # carrier trial must have exactly one V48 operator trial.
    comparable_trials=[r for r in trial_rows if int(numeric(r,'iteration'))>1]
    if len(attempt_rows)!=len(comparable_trials):
        raise ValueError("V48 path not evaluated for every comparable V47 carrier trial")
    for row,v47 in zip(attempt_rows,comparable_trials):
        for key in ('iteration','retry','beta'):
            if numeric(row,key)!=numeric(v47,key):
                raise ValueError("V47 and V48 attempt identity mismatch")
        if numeric(row,'previous_identity_l1')>1e-10 or numeric(row,'previous_identity_linf')<0:
            raise ValueError("previous committed RC flux identity broken")
        if row.get('effective_rAU_rederived')!='1' or row.get('exact_physical_residual_oracle')!='0':
            raise ValueError("V48 operator provenance missing")
    partial={(int(numeric(row,'iteration')),int(numeric(row,'retry')))
             for row in attempt_rows if numeric(row,'beta')<1.0}
    native={(int(numeric(row,'iteration')),int(numeric(row,'retry')))
             for row in native_rows}
    if native != partial:
        raise ValueError("missing, extra or duplicated native coefficient discrepancy")
    for row in native_rows:
        for key in ('relative_l1','absolute_linf'):
            if numeric(row,key)<0:
                raise ValueError('negative native coefficient discrepancy')
        if row.get('native_operator_qualified')!='0':
            raise ValueError('homotopy native operator promoted without oracle')
    commits=[]
    last_merit=math.inf
    for row in commit_rows:
        iteration=int(numeric(row,'iteration'))
        merit=numeric(row,'merit')
        if merit<0 or (commits and iteration<=commits[-1]['iteration']):
            raise ValueError('invalid committed iteration or residual')
        comps=[numeric(row,key) for key in ('U','p_continuity','k','omega','nut_closure')]
        if any(c<0 for c in comps) or abs(merit-max(comps))>1e-8*max(1.,merit):
            raise ValueError('global physical components inconsistent with merit')
        if commits and not merit<last_merit:
            raise ValueError('V14 monotonicity violated by committed state')
        if row.get('qualified')!='0':
            raise ValueError('experimental homotopy cannot claim physical qualification')
        if iteration>1 and not any(int(numeric(a,'iteration'))==iteration for a in attempt_rows):
            raise ValueError('committed state has no V48 operator trial (except first)')
        commits.append({'iteration':iteration,'merit':merit,
                        'converged_flag':row.get('converged','?'),
                        'continuity':comps[1]})
        last_merit=merit
    if len(commits)<min_commits:
        raise ValueError(f'only {len(commits)} committed states (<{min_commits})')
    aborted='CFDX production solver error:' in log or bool(re.search(r'nonlinear iteration \d+ rejected',log)) and 'Converged NO' not in log
    finished=bool(re.search(r'Converged (YES|NO) Iterations \d+',log))
    status=('ABORTED_UNQUALIFIED' if aborted else
            'INCOMPLETE_UNQUALIFIED' if finished else 'TRUNCATED_UNQUALIFIED')
    return {'status':status,'physical_qualification':False,
            'path_dependent_rau':True,'commits':len(commits),
            'last':commits[-1],'native_gap_samples':len(native_rows),
            'v48_coefficient_identity_trials':len(attempt_rows),
            'attempts':len(trial_rows),'partial_beta_attempts':len(partial),
            'default_behavior_modified':False}

def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log',required=True,type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--min-commits',type=int,default=1)
    args=p.parse_args()
    result=inspect(args.log.read_text(errors='replace'),args.min_commits)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))
if __name__=='__main__':
    main()
