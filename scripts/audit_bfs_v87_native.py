"""Fail-closed V87 pre-SST and post-SST audit. Not a convergence certificate."""
import argparse, json, re
from pathlib import Path

def parse_line(line):
    return dict(re.findall(r'([A-Za-z_]+)=([^\s]+)',line))
def audit(text):
    pending=None
    passed=0
    rejected=0
    retries=0
    committed=[]
    v86_checks=0
    events=[]
    for lineno,line in enumerate(text.splitlines(),1):
        if line.startswith('BFS_EARLY_FIELD_TRIAL_V87 '):
            if pending is not None:
                raise ValueError(f'Previous early trial never completed before line {lineno}')
            d=parse_line(line)
            if d.get('qualified')!='0' or d.get('continue_sst') not in ('0','1'):
                raise ValueError('Invalid or falsely certified V87 trial')
            ratio=float(d['max_speed_over_inlet'])
            if not (0<=ratio<1e6):raise ValueError('Invalid velocity ratio')
            if d['continue_sst']=='0' and (ratio<=2.0 or d.get('reason')!='velocity_envelope'):
                raise ValueError('Incoherent rejected velocity trial')
            if d['continue_sst']=='1' and ratio>2.0:
                raise ValueError('Unsafe pre-SST velocity accepted')
            pending=('rejected' if d['continue_sst']=='0' else 'passed',d['iteration'])
            events.append({'line':lineno,'iteration':int(d['iteration']),
                           'result':pending[0],'speed_ratio':ratio})
        elif 'phase=SST_TRANSPORT_BEGIN' in line and line.startswith('BFS_PHASE_V10 '):
            d=parse_line(line)
            if pending is not None:
                if pending[1]!=d['iteration'] or pending[0]!='passed':
                    raise ValueError('SST ran after rejected V87 preflight')
                passed+=1
                pending=None
        elif line.startswith('CFDX nonlinear iteration ') and 'rejected' in line:
            n=re.search(r'iteration (\d+)',line)
            if pending is not None:
                if pending[0]!='rejected' or not n or pending[1]!=n.group(1):
                    raise ValueError('Missing immediate nonlinear rollback after preflight rejection')
                rejected+=1
                pending=None
            retries+=1
        elif line.startswith('BFS_FIELD_TRUST_TRIAL_V86 '):
            d=parse_line(line)
            if d.get('qualified')!='0':raise ValueError('False V86 physical certificate')
            v86_checks+=1
        elif line.startswith('GLOBAL_PROGRESS_COMMITTED_V85 '):
            d=parse_line(line)
            if d.get('qualified')!='0' or d.get('converged')!='0':
                raise ValueError('False final certificate in provisional V85 progress')
            committed.append(int(d['iteration']))
    # Timed-out logs may end after preflight on an unfinished SST solve;
    # this may not be falsely labelled a pass. Record as incomplete.
    if any(b<=a for a,b in zip(committed,committed[1:])):
        raise ValueError('Non-monotonic or duplicate committed iterations')
    return {'schema':'CFDX-V87-PRESST-EVIDENCE/1','early_trials':len(events),
            'early_rejections_with_rollback':rejected,'early_passes_to_sst':passed,
            'incomplete_last_preflight':pending is not None,
            'nonlinear_retries':retries,'post_sst_v86_checks':v86_checks,
            'committed_iterations':committed,'all_committed_provisional':True,
            'physical_convergence_certified':False,'events':events}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('native_stderr',type=Path)
    p.add_argument('--output',type=Path)
    a=p.parse_args();r=audit(a.native_stderr.read_text(errors='replace'))
    print(json.dumps(r,indent=2))
    if a.output:a.output.write_text(json.dumps(r,indent=2))