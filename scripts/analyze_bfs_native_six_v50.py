"""Fail-closed forensic verification of V50's *diagnostic* six-field vector.

All output is physical qualification = false, even for a closed native face
flux; a two-way equality of scalar and vector observables is NOT convergence.
"""
from __future__ import annotations
import argparse
import json
import math
import pathlib
import re
from typing import Any

TOKENS = re.compile(r'([^\s=]+)=([^\s]+)')

def parse_line(line: str, kind: str) -> dict[str, str]:
    if not line.startswith(kind + ' '):
        raise ValueError('incorrect forensic record type')
    result = dict(TOKENS.findall(line[len(kind):]))
    if not result:
        raise ValueError('empty forensic record')
    return result

def integer(record: dict[str, str], key: str) -> int:
    try:
        value=int(record[key])
    except (ValueError, KeyError) as exc:
        raise ValueError(f'missing/invalid {key}') from exc
    return value

def finite(record: dict[str, str], key: str) -> float:
    try:
        value=float(record[key])
    except (ValueError, KeyError) as exc:
        raise ValueError(f'missing/invalid {key}') from exc
    if not math.isfinite(value):
        raise ValueError(f'nonfinite {key}')
    return value

def verify(text: str, *, min_closed: int=1, min_commits: int=1) -> dict[str, Any]:
    attempted: dict[tuple[int,int],dict[str,str]]={}
    native: dict[tuple[int,int],dict[str,str]]={}
    commits: dict[int,dict[str,str]]={}
    v48_commits: dict[int,dict[str,str]]={}
    last_trial: dict[int,tuple[int,int]]={}
    final_gates: dict[int,dict[str,str]]={}
    for line in text.splitlines():
        if line.startswith('NATIVE_SIX_VECTOR_V50 '):
            r=parse_line(line,'NATIVE_SIX_VECTOR_V50')
            key=(integer(r,'iteration'),integer(r,'retry'))
            if key in attempted or key[0]<1 or key[1]<0:
                raise ValueError('duplicate or invalid V50 trial')
            attempted[key]=r
            last_trial[key[0]]=key
        elif line.startswith('NATIVE_RC_ORACLE_V49 '):
            r=parse_line(line,'NATIVE_RC_ORACLE_V49')
            key=(integer(r,'iteration'),integer(r,'retry'))
            if key in native:
                raise ValueError('duplicate V49 trial')
            native[key]=r
        elif line.startswith('NATIVE_SIX_COMMITTED_V50 '):
            r=parse_line(line,'NATIVE_SIX_COMMITTED_V50')
            it=integer(r,'iteration')
            if it in commits:
                raise ValueError('duplicate V50 commit')
            commits[it]=r
        elif line.startswith('NATIVE_SIX_FINAL_GATE_V50 '):
            r=parse_line(line,'NATIVE_SIX_FINAL_GATE_V50')
            it=integer(r,'iteration')
            if it in final_gates:
                raise ValueError('duplicate native final-stop gate')
            final_gates[it]=r
        elif line.startswith('GLOBAL_NONLINEAR_COMMITTED_V13 '):
            r=parse_line(line,'GLOBAL_NONLINEAR_COMMITTED_V13')
            it=integer(r,'iteration')
            if it in v48_commits:
                raise ValueError('duplicate original commit')
            v48_commits[it]=r
    if not attempted or not v48_commits or not commits:
        raise ValueError('missing trial/native/committed evidence')
    if set(commits)!=set(v48_commits):
        raise ValueError('V50 commits do not match the production transactions')
    if set(attempted)!=set(native):
        raise ValueError('V50 trial unmatched to V49 native fixed point')
    closed=0
    for key, v50 in attempted.items():
        v49=native[key]
        if integer(v50,'qualified') or integer(v50,'committed') or \
           integer(v49,'qualified') or integer(v49,'physical_converged'):
            raise ValueError('false physical qualification or tentative commit')
        is_closed=integer(v50,'closed')
        if is_closed not in (0,1) or is_closed!=integer(v49,'closed'):
            raise ValueError('contradictory V49/V50 closure')
        if not is_closed:
            continue
        closed+=1
        cells=integer(v50,'cells')
        if cells<=0 or integer(v50,'entries')!=6*cells:
            raise ValueError('incomplete six-field vector')
        if integer(v50,'same_native_flux')!=1 or integer(v50,'parity')!=1 or \
           integer(v50,'experimental')!=1:
            raise ValueError('missing common operator or parity evidence')
        fingerprint=integer(v50,'fingerprint')
        if fingerprint<0 or fingerprint>=2**64:
            raise ValueError('invalid vector fingerprint')
        for label in ('U','continuity','k','omega'):
            a=finite(v50,label)
            b=finite(v49,label)
            if a<0 or b<0 or abs(a-b)>1e-9*max(1.0,a,b):
                raise ValueError(f'scalar/vector {label} mismatch')
        merit=finite(v50,'merit')
        if merit<0 or abs(merit-max(finite(v50,k) for k in
                                  ('U','continuity','k','omega')))>1e-9*max(1.0,merit):
            raise ValueError('invalid six-field merit')
    if closed<min_closed or len(commits)<min_commits:
        raise ValueError('insufficient closed six-field evidence')
    rows=[]
    for it in sorted(commits):
        c=commits[it]
        if integer(v48_commits[it],'converged')!=0 or \
           integer(v48_commits[it],'qualified')!=0:
            raise ValueError('production log falsely claims global physical convergence')
        key=last_trial.get(it)
        if key is None or integer(attempted[key],'closed')!=1:
            raise ValueError('committed state does not have a closed native vector')
        trial=attempted[key]
        if integer(c,'committed')!=1 or integer(c,'qualified')!=0 or \
           integer(c,'physical_converged')!=0:
            raise ValueError('invalid V50 commit certificate')
        for field in ('fingerprint','cells','entries'):
            if integer(c,field)!=integer(trial,field):
                raise ValueError('committed fingerprint or shape differs from trial')
        for field in ('merit','U','continuity','k','omega'):
            if abs(finite(c,field)-finite(trial,field))>1e-9*max(
                    1.0,abs(finite(c,field)),abs(finite(trial,field))):
                raise ValueError('committed value differs from tentative state')
        if it in final_gates:
            gate=final_gates[it]
            allowed=integer(gate,'eligible')
            if integer(gate,'experimental')!=1 or \
               integer(gate,'qualified')!=0 or allowed not in (0,1):
                raise ValueError('invalid native final-stop gate')
            tolerance=finite(gate,'tolerance')
            gap=finite(gate,'gap')
            if tolerance<=0 or gap<0 or gap!=finite(c,'flux_gap_l1'):
                raise ValueError('native final-stop flux gap mismatch')
            expected=(finite(c,'merit')<=tolerance and gap<=tolerance)
            # The extra dormant-momentum global check may only turn an
            # apparently eligible scalar merit into a rejection.
            if allowed and not expected:
                raise ValueError('false native final-stop acceptance')
            if integer(gate,'native_available')!=1:
                raise ValueError('native stop has no complete vector')
        rows.append({'iteration':it,'native_merit':finite(c,'merit'),
                     'native_k':finite(c,'k'),'native_omega':finite(c,'omega'),
                     'production_merit':finite(v48_commits[it],'merit'),
                     'production_k':finite(v48_commits[it],'k'),
                     'production_omega':finite(v48_commits[it],'omega'),
                     'fingerprint':integer(c,'fingerprint')})
    if final_gates and set(final_gates)!=set(commits):
        raise ValueError('incomplete native final-stop instrumentation')
    return {'status':'DIAGNOSTIC_ONLY_PHYSICAL_UNQUALIFIED','qualified':False,
            'closed_trials':closed,'all_trials':len(attempted),
            'committed_outer_states':len(commits),'cells':integer(commits[min(commits)],'cells'),
            'states':rows}

def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--log',type=pathlib.Path,required=True)
    ap.add_argument('--output',type=pathlib.Path,required=True)
    ap.add_argument('--min-closed',type=int,default=1)
    ap.add_argument('--min-commits',type=int,default=1)
    args=ap.parse_args()
    try:
        result=verify(args.log.read_text(errors='replace'),
                      min_closed=args.min_closed,min_commits=args.min_commits)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        print(json.dumps(result,sort_keys=True))
    except (OSError,ValueError) as exc:
        print('V50 rejected evidence:',exc)
        return 2
    return 0
if __name__=='__main__':
    raise SystemExit(main())
