@@ -0,0 +1,167 @@
#!/usr/bin/env python3
"""Fail-closed read-only BFS152 V66 evidence auditor.

Consumes ONLY genuine C++ V65/V66 solver trace lines. The result is a
forensic report, never a physical or NASA validation certificate. A rejected
V65 step is still rejected even when some V66 consistency checks succeed.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re

KV = re.compile(r'([A-Za-z][A-Za-z_0-9]*)=([^\s]+)')
MARKERS = {
    'v65': 'GLOBAL_SAME_OPERATOR_TRIAL_V65 ',
    'carrier': 'BFS_CARRIER_COHERENCE_V66 ',
    'omega': 'BFS_OMEGA_HOTSPOT_V66 ',
    'budget': 'BFS_OMEGA_BUDGET_V66 ',
}

def read_rows(text, prefix):
    return [dict(KV.findall(line)) for line in text.splitlines()
            if line.startswith(prefix)]

def number(row,key,nonnegative=False):
    if key not in row:raise ValueError(f'missing {key}')
    v=float(row[key])
    if not math.isfinite(v) or (nonnegative and v<0):
        raise ValueError(f'invalid {key}')
    return v

def integer(row,key):
    if key not in row:raise ValueError(f'missing {key}')
    v=int(row[key]);
    if v<0:raise ValueError(f'invalid {key}')
    return v

def trial_key(row):return integer(row,'iteration'),integer(row,'retry')

def near(a,b,rtol=1e-9,atol=1e-12):return abs(a-b)<=atol+rtol*max(abs(a),abs(b))

def audit(text):
    groups={k:read_rows(text,v) for k,v in MARKERS.items()}
    v65={trial_key(row):row for row in groups['v65']}
    if len(v65)!=len(groups['v65']):raise ValueError('duplicate V65 trace')
    carriers={}
    for row in groups['carrier']:
        key=trial_key(row)
        if key in carriers:raise ValueError('duplicate V66 carrier trace')
        carriers[key]=row
    expected={key for key in v65 if key[0] in (22,23)}
    if not expected or not any(x[0]==23 for x in expected):
        raise ValueError('V66 window not reached (iteration 23 required)')
    if expected!=set(carriers):
        raise ValueError('missing/unmatched V66 or V65 trial evidence')
    hotspots={}
    for row in groups['omega']:
        key=trial_key(row)
        hotspots.setdefault(key,[]).append(row)
    if set(hotspots)!=expected:raise ValueError('missing/alien V66 omega hotspot')
    budgets={}
    for row in groups['budget']:
        key=trial_key(row)
        if key in budgets:raise ValueError('duplicate V66 omega budget')
        budgets[key]=row
    probes={k for k in expected if k[1]==0 or k[1]>=8}
    if set(budgets)!=probes:raise ValueError('missing/alien V66 row budget')
    result=[]
    for key in sorted(expected):
        v=v65[key];c=carriers[key];h=hotspots[key]
        if len(h)!=3:raise ValueError('must report three omega rankings per trial')
        if not near(number(v,'beta'),number(c,'beta')):
            raise ValueError('V65/V66 beta mismatch')
        for metric in ('old_carrier_gap','trial_carrier_gap','rau_change',
                       'nut_change','old_rau_carrier_gap','trial_rau_carrier_gap',
                       'old_nut_carrier_gap','trial_nut_carrier_gap','old_carrier_div_gap','trial_carrier_div_gap'):
            number(c,metric,True)
        if integer(c,'old_rau_carrier_available') not in (0,1):
            raise ValueError('invalid production rAU evidence')
        if not near(number(c,'old_carrier_gap'),number(v,'old_authoritative_phi_gap'),atol=1e-12) or \
           not near(number(c,'trial_carrier_gap'),number(v,'trial_authoritative_phi_gap'),atol=1e-12):
            raise ValueError('V65/V66 carrier gap mismatch')
        omega=[]
        for row in h:
            if integer(row,'cell')>=152:raise ValueError('BFS152 invalid cell id')
            for k in ('x','y','z','old_raw','trial_raw','delta',
                      'omega_old','omega_trial','k_old','k_trial','nut_old','nut_trial'):
                number(row,k)
            if not near(number(row,'trial_raw')-number(row,'old_raw'),number(row,'delta')):
                raise ValueError('fake omega signed-row delta')
            omega.append({'cell':int(row['cell']),
                          'raw_old':float(row['old_raw']),
                          'raw_trial':float(row['trial_raw']),
                          'delta':float(row['delta']),
                          'position':[float(row[t]) for t in ('x','y','z')]})
        if omega[0]['cell'] != integer(c,'omega_change_cell') or \
           omega[1]['cell'] != integer(c,'omega_trial_cell') or \
           omega[2]['cell'] != integer(c,'omega_old_cell'):
            raise ValueError('hotspot ranks inconsistent with carrier trace')
        row_budget=None
        if key in budgets:
            b=budgets[key]
            if integer(b,'cell') != omega[1]['cell']:
                raise ValueError('SST budget not on worst trial omega row')
            dterms=['convection','diffusion','production','destruction','cross']
            terms={term:number(b,term) for term in dterms}
            total=sum(terms.values())
            if not near(total,number(b,'sum'),atol=1e-9) or \
               not near(total,number(b,'residual_delta'),atol=1e-9) or \
               not near(number(b,'trial_row')-number(b,'old_row'),number(b,'residual_delta'),atol=1e-9):
                raise ValueError('SST term decomposition does not close')
            for field in ('closure','oracle_gap_old','oracle_gap_trial'):
                if abs(number(b,field))>1e-8*max(1.,abs(number(b,'trial_row')),abs(number(b,'old_row'))):
                    raise ValueError('unclosed SST FV row')
            number(b,'old_wall_flux');number(b,'trial_wall_flux')
            if integer(b,'owner_constraint') not in (0,1):
                raise ValueError('invalid SST owner constraint')
            row_budget=terms
        if integer(c,'committed')!=0 or integer(c,'physical_qualified')!=0 or \
           integer(v,'committed')!=0 or integer(v,'qualified')!=0:
            raise ValueError('forensic probe promoted to physical PASS')
        result.append({'iteration':key[0],'retry':key[1],
                       'beta':number(c,'beta'),
                       'carrier_gap':number(c,'trial_carrier_gap'),
                       'carrier_face':integer(c,'worst_face'),
                       'rAU_relative_change':number(c,'rau_change'),
                       'nut_relative_change':number(c,'nut_change'),
                       'omega':omega,'omega_budget_delta':row_budget,
                       'v65_accepted':bool(integer(v,'accepted'))})
    return {'schema':'CFDX-BFS-V66-CAUSAL/1','cases':result,
            'trial_count':len(result),'budget_count':len(budgets),
            'physical_qualification':False,'convergence_claim':False}

def check_unmodified_v65_trace(current,reference):
    """Bit-for-bit evidence that V66 did not modify accepted/rejected V65 updates."""
    prefixes=('GLOBAL_SAME_OPERATOR_TRIAL_V65 ',
              'GLOBAL_SAME_OPERATOR_COMMITTED_V65 ',
              'GLOBAL_NONLINEAR_COMMITTED_V13 ')
    comparison={}
    for prefix in prefixes:
        a=[line for line in current.splitlines() if line.startswith(prefix)]
        b=[line for line in reference.splitlines() if line.startswith(prefix)]
        if not a or a!=b:
            raise ValueError('V66 changes baseline V65 physical trace: '+prefix.strip())
        comparison[prefix.strip()]=len(a)
    return comparison

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--reference-log',type=Path,help='optional original V65 trace for exact parity')
    a=p.parse_args()
    raw=a.log.read_text()
    result=audit(raw)
    if a.reference_log:
        result['reference_v65_parity']=check_unmodified_v65_trace(
            raw,a.reference_log.read_text())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('V66_CAUSAL_AUDIT_OK trials=',result['trial_count'],
          'budgets=',result['budget_count'],'physical_qualification=0')
if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as ex:raise SystemExit('V66_AUDIT_FAIL: '+str(ex))
