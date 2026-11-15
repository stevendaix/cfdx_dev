"""Fail-closed read-only BFS152 audit of sole-phi SST FV substitution.

This is NOT a convergence certificate. Evidence only passes if the same-state,
convection-only identities close and the V65 decisions remain identical.
"""
import argparse
import json
import math
from pathlib import Path

PREFIX = 'BFS_SST_FLUX_SWAP_V68 '
COMPARE_PREFIXES = (
    'GLOBAL_SAME_OPERATOR_TRIAL_V65 ',
    'CFDX nonlinear iteration ',
)
REQUIRED = (
    'iteration', 'retry', 'face', 'phi_carrier', 'phi_physical',
    'face_omega_owner', 'face_omega_neighbour', 'omega_phys_linf',
    'omega_carrier_linf', 'omega_effect_fraction', 'omega_delta_linf',
    'omega_worst_phys', 'omega_worst_carrier', 'omega_worst_delta',
    'omega_phys_cell0', 'omega_carrier_cell0', 'omega_phys_cell12',
    'omega_carrier_cell12', 'k_phys_linf', 'k_carrier_linf',
    'k_delta_linf', 'face_identity_gap', 'source_diffusion_gap',
    'same_six_state', 'sole_phi_changed', 'committed', 'qualified',
)

def read_records(text):
    rows=[]
    for line in text.splitlines():
        if not line.startswith(PREFIX):
            continue
        attrs={}
        for part in line[len(PREFIX):].split():
            if '=' not in part:
                raise ValueError('malformed V68 evidence record')
            key, value=part.split('=',1)
            if key in attrs:
                raise ValueError('duplicate key '+key)
            attrs[key]=value
        if any(key not in attrs for key in REQUIRED):
            raise ValueError('incomplete flux swap record')
        numeric={key: float(value) for key, value in attrs.items()}
        if not all(math.isfinite(v) for v in numeric.values()):
            raise ValueError('nonfinite flux swap evidence')
        if (int(numeric['iteration'])!=23 or
            int(numeric['retry']) not in (0,8,10) or
            int(numeric['face'])!=3):
            raise ValueError('off-scope V68 evidence')
        if (numeric['same_six_state']!=1 or numeric['sole_phi_changed']!=1 or
            numeric['committed']!=0 or numeric['qualified']!=0):
            raise ValueError('flux substitution is not read-only / same-state')
        if (numeric['face_identity_gap']>1.e-9 or
            numeric['source_diffusion_gap']>1.e-8 or
            numeric['face_identity_gap']<0 or
            numeric['source_diffusion_gap']<0):
            raise ValueError('V68 independent FV identity does not close')
        if any(numeric[name]<0 for name in ('omega_phys_linf','omega_carrier_linf',
                                           'omega_delta_linf','k_phys_linf','k_carrier_linf')):
            raise ValueError('negative absolute raw FV residual')
        expected=(numeric['omega_phys_linf']-numeric['omega_carrier_linf'])/max(numeric['omega_phys_linf'],1.e-300)
        if abs(expected-numeric['omega_effect_fraction'])>1e-8:
            raise ValueError('incorrect omega change fraction')
        rows.append(numeric)
    if len(rows)!=3 or sorted(int(v['retry']) for v in rows)!=[0,8,10]:
        raise ValueError('missing/duplicated V68 full/small-step probes')
    return sorted(rows,key=lambda v:v['retry'])

def v65_decisions(text):
    return [line for line in text.splitlines() if line.startswith(COMPARE_PREFIXES)]

def audit(experimental, baseline):
    rows=read_records(experimental)
    decisions=v65_decisions(experimental)
    ref=v65_decisions(baseline)
    if not ref or decisions != ref:
        raise ValueError('V68 changed/missed V65 line-search decisions')
    accepted=sum('accepted=1' in line for line in decisions
                 if line.startswith('GLOBAL_SAME_OPERATOR_TRIAL_V65 '))
    if accepted!=22 or sum(line.startswith('GLOBAL_SAME_OPERATOR_TRIAL_V65 ')
                           for line in decisions)!=149:
        raise ValueError('V65 reference campaign not fully reproduced')
    if 'CFDX production solver error: V65 state-only physical evaluation failed' not in experimental:
        raise ValueError('V68 must not certify nonconverged physical campaign')
    return {'schema':'cfdx-v68-same-six-state-sole-flux-swap-1',
            'read_only':True,'physical_qualification':False,
            'v65_replay_equal':True,'v65_trials':149,'v65_accepted':22,
            'v65_failure_iteration':23,'flux_swap_probes':rows,
            'interpretation':'omega_effect_fraction is the relative change in raw infinity norm after swapping ONLY mass flux; it is not a global merit improvement or a converged SST solution'}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--log',type=Path,required=True)
    ap.add_argument('--baseline',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    report=audit(a.log.read_text(),a.baseline.read_text())
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('V68_SST_FLUX_SWAP_EVIDENCE_PASS')
if __name__=='__main__':main()