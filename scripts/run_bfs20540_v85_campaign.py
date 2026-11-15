#!/usr/bin/env python3
"""Run bounded V85 native BFS campaign with fail-closed evidence summary.

Native progress does NOT imply convergence or a valid NASA comparison.
The operator's maximum outer iteration count is also capped in the case H5.
"""
import argparse,csv,json,os,re,subprocess,time,hashlib
from pathlib import Path

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def records(log,prefix):
    if not Path(log).exists():return []
    rows=[]
    for line in Path(log).read_text(errors='replace').splitlines():
        if line.startswith(prefix+' '):
            rows.append(dict(re.findall(r'([a-zA-Z_]+)=([^\s]+)',line)))
    return rows

def run(binary,case,out,max_seconds=900,limit=100):
    binary,case,out=Path(binary),Path(case),Path(out)
    if not binary.is_file() or not case.is_file():raise FileNotFoundError('Missing binary or case')
    if not 1<=limit<=100:raise ValueError('Upper iteration bound must be <=100')
    out.mkdir(parents=True,exist_ok=True)
    cmd=[str(binary.resolve()),'--mesh',str(case.resolve()),'--experimental-sst',
         '--iterations',str(limit),'--output-dir',str(out.resolve())]
    stdout=out/'native.stdout';stderr=out/'native.stderr';start=time.monotonic()
    env=dict(os.environ,OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS','2'),
             CFDX_V80_NATIVE_FACE_DEFECT='1',CFDX_SST_DIAGNOSTICS='1',CFDX_BFS_STAGE_TRACE='1')
    with stdout.open('w') as so,stderr.open('w') as se:
        proc=subprocess.Popen(cmd,stdout=so,stderr=se,env=env)
        try:code=proc.wait(timeout=max_seconds);timed_out=False
        except subprocess.TimeoutExpired:
            proc.terminate()
            try:proc.wait(timeout=5)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
            code=124;timed_out=True
    progress=records(stderr,'GLOBAL_PROGRESS_COMMITTED_V85')
    candidates=records(stderr,'GLOBAL_PROGRESS_TRIAL_V85')
    physical=records(stderr,'BFS_FIELD_COMMITTED_V85')
    iterations=[int(x['iteration']) for x in progress]
    if any(b<=a for a,b in zip(iterations,iterations[1:])) or any(i>limit for i in iterations):
        raise RuntimeError('Corrupted V85 committed iteration sequence')
    pmap={int(x['iteration']):x for x in physical}
    if set(iterations)!=set(pmap):raise RuntimeError('Missing per-state field proof')
    if any(r.get('qualified')!='0' for r in progress):raise RuntimeError('False V85 qualification marker')
    vtu=sorted(out.glob('result_*.vtu'))
    for file in vtu:
        it=int(file.stem.rsplit('_',1)[-1])
        if it not in iterations or (it>10 and it%10!=0 and it!=limit):
            raise RuntimeError('VTU not associated with a permitted native committed state')
    columns=['iteration','U','U_global','p_continuity','k','omega','merit',
             'reverse_ux_cells','ux_min','ux_max','pressure_min_pa','pressure_max_pa','pressure_span_pa']
    with (out/'accepted_100_monitor.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,columns);writer.writeheader()
        for r in progress:
            row={key:r.get(key,'') for key in columns}
            for key in ('reverse_ux_cells','ux_min','ux_max','pressure_min_pa','pressure_max_pa','pressure_span_pa'):
                row[key]=pmap[int(r['iteration'])].get(key,'')
            writer.writerow(row)
    report={'schema':'CFDX_V85_BFS20540_NATIVE_CAMPAIGN_1','status':'not_finally_qualified',
        'max_outer_requested':limit,'outer_states_accepted':len(progress),
        'last_accepted_iteration':max(iterations,default=0),
        'native_process_code':code,'timed_out':timed_out,
        'elapsed_seconds':round(time.monotonic()-start,2),'command':cmd,
        'committed_records':len(progress),'rejected_proposals':sum(x.get('accepted')=='0' for x in candidates),
        'native_vtu_files':[{'name':f.name,'sha256':sha(f)} for f in vtu],
        'native_stderr_sha256':sha(stderr),'input_case_sha256':sha(case),
        'targets':[10,20,50,100],
        'physical_validation':'none: provisional SST, final V8 criteria must be independently satisfied'}
    (out/'campaign_report.json').write_text(json.dumps(report,indent=2))
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--case',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--max-runtime',type=int,default=900)
    p.add_argument('--max-outer',type=int,default=100)
    a=p.parse_args();print(json.dumps(run(a.binary,a.case,a.out,a.max_runtime,a.max_outer),indent=2))