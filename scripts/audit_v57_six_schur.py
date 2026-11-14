"""Fail-closed CFDX V57: six native Schur methods on one 4N + 6N state.

This is a *negative qualification* reporter. Even algebraically clean Krylov
snapshots are NOT proof of physical convergence or OpenFOAM parity. Matrix
fingerprints reported by CFDX V40 are FNV, **not** independent SHA-256.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

METHODS = ('bfbt', 'block_local', 'lsc', 'simple', 'simplec', 'pcd')
LINEAR_LIMIT = 0.3  # unchanged V52 experimental target
KRYLOV_LIMIT = 1e-6  # additional V57 consistency check, never a physical gate

def tokens(line):
    return dict(re.findall(r'([a-zA-Z_][a-zA-Z_0-9]*)=([^\s]+)', line))

def as_float(v, name):
    try:
        out = float(v)
    except (ValueError, TypeError):
        raise ValueError(f'{name}: invalid numeric field') from None
    if not math.isfinite(out):
        raise ValueError(f'{name}: nonfinite field')
    return out

def as_int(v, name):
    try:
        out=int(v)
    except (ValueError, TypeError):
        raise ValueError(f'{name}: invalid integer field') from None
    if str(out)!=str(v):
        raise ValueError(f'{name}: noncanonical integer field')
    return out

def analyze_log(method: str, path: Path):
    if method not in METHODS:
        raise ValueError('unknown native Schur method')
    raw=path.read_bytes()
    lines=raw.decode('utf-8').splitlines()
    selected=lambda prefix: [tokens(l) for l in lines if l.startswith(prefix+' ')]
    matr=selected('COUPLED_MATRIX_FINGERPRINT_V40')
    if len(matr)!=1 or matr[0].get('fnv64_not_sha256')!='1':
        raise ValueError(f'{method}: missing authoritative 4N fingerprint')
    sig=tuple(matr[0].get(k) for k in ('A','b','x0','rows','nnz'))
    if None in sig or sig[3]!='608' or as_int(sig[4],'nnz')<=0:
        raise ValueError(f'{method}: incompatible 4N matrix')
    summary=selected('FPHYS_JFNK_V52')
    snapshots=selected('FPHYS_KRYLOV_TRUTH_V57')
    if method=='pcd':
        if summary or snapshots or not selected('COUPLED_LINEAR_UNCHANGED_RETRY_V37') or \
           not any('coupled momentum-continuity solve did not converge' in l for l in lines):
            raise ValueError('PCD result must expose its 4N failure, not a fallback')
        return {'method':method,'status':'BLOCKED_PCD_4N','matrix':sig,
                'log_sha256':hashlib.sha256(raw).hexdigest(), 'krylov_directions':0}
    if len(summary)!=1 or not snapshots:
        raise ValueError(f'{method}: missing 6N true-residual audit')
    result=summary[0]
    if result.get('selected_schur')!=method or result.get('frozen_prestep')!='1' or \
       result.get('v53_global_block')!='1' or result.get('central_jv')!='1' or \
       result.get('private_accepted')!='0' or result.get('production_committed')!='0' or \
       result.get('qualified')!='0' or result.get('v57_hybrid_jv_calls') is None:
        raise ValueError(f'{method}: solver/state/op/operator qualification contract mismatch')
    if as_int(result.get('v53_transport_nnz'),'nnz') != 1720:
        raise ValueError(f'{method}: SST spatial block not identical')
    if len(snapshots)!=as_int(result.get('v57_true_krylov_checks'),'directions'):
        raise ValueError(f'{method}: missing direction certificates')
    if len(snapshots)!=16:
        raise ValueError(f'{method}: Krylov budgets are not common (16)')
    checks=[]
    for i, row in enumerate(snapshots,1):
        if as_int(row.get('direction'),'direction')!=i or row.get('production_committed')!='0':
            raise ValueError(f'{method}: out-of-order or committed Krylov direction')
        if as_float(row.get('consistency_limit'),'consistency_limit')!=KRYLOV_LIMIT:
            raise ValueError(f'{method}: mismatched Krylov guard')
        values={k:as_float(row.get(k),k) for k in ('arnoldi','true','gap','nonadditivity')}
        if any(v<0 for v in values.values()):
            raise ValueError(f'{method}: negative Krylov defect')
        if abs(values['true']-values['arnoldi'])>values['gap']+1e-13:
            raise ValueError(f'{method}: inconsistent Arnoldi gap')
        certified=(values['gap']<=KRYLOV_LIMIT and
                   values['nonadditivity']<=KRYLOV_LIMIT)
        if row.get('qualified')!=str(int(certified)):
            raise ValueError(f'{method}: false intermediate direction certification')
        checks.append(values)
    true_final=as_float(result.get('linear_true_relative'),'linear_true_relative')
    if abs(true_final-checks[-1]['true'])>1e-10:
        raise ValueError(f'{method}: true residual final not reproduced')
    if result.get('v57_hybrid_consistent')!='0' or result.get('v57_krylov_consistent')!='0':
        raise ValueError(f'{method}: unqualified direction promoted')
    conv_gap=as_float(result.get('v57_hybrid_convection_gap'),'conv_gap')
    nonsmooth=as_int(result.get('v57_hybrid_nonsmooth_faces'),'nonsmooth')
    if nonsmooth<=0 and conv_gap<=1e-3:
        raise ValueError(f'{method}: expected branch proof not present')
    if true_final<=LINEAR_LIMIT:
        raise ValueError(f'{method}: unexpected linear success requires fresh qualification')
    if not any(l.strip()=='Converged NO Iterations 1' for l in lines):
        raise ValueError(f'{method}: unqualified physical termination missing')
    return {'method':method,'status':'BLOCKED_NONSMOOTH_CONVECTION',
            'matrix':sig,'common_state':result.get('common_state_fingerprint'),
            'prestep_merit':as_float(result.get('merit_before'),'merit_before'),
            'true_linear_relative':true_final,
            'max_arnoldi_gap':max(c['gap'] for c in checks),
            'max_nonadditivity':max(c['nonadditivity'] for c in checks),
            'max_convective_gap':conv_gap,'nonsmooth_component_face_events':nonsmooth,
            'krylov_directions':len(checks),
            'log_sha256':hashlib.sha256(raw).hexdigest()}

def audit_six(folder:Path):
    evidence=[analyze_log(m,folder/f'six_schur_{m}.log') for m in METHODS]
    if len({tuple(row['matrix']) for row in evidence})!=1:
        raise ValueError('six Schur used different frozen 4N matrices')
    usable=evidence[:-1]
    if any(not row['common_state'] for row in usable) or \
       len({row['common_state'] for row in usable})!=1 or \
       len({row['prestep_merit'] for row in usable})!=1:
        raise ValueError('six Schur do not share the same frozen 6N state')
    return {'version':'CFDX_V57','qualified':False,
            'status':'BLOCKED_GLOBAL_JV_AND_PCD',
            'same_4n_FNV_matrix':True, 'same_6n_initial_state':True,
            'note':'FNV64 matrix check; physical nonlinear convergence not certified',
            'common_state_fingerprint':usable[0]['common_state'],
            'methods':evidence}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--logs-dir',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    report=audit_six(args.logs_dir)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(f"V57_SIX_SCHUR {report['status']} physical_qualified=0")
if __name__=='__main__':main()
