#!/usr/bin/env python3
"""Run ONE bounded, private first SST transport on a given CFDX HDF5 case.

Stops immediately after the first SST_TRANSPORT_END (even if CFDX intends
nonlinear retries). It is a DIAGNOSTIC, never a convergence certification.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1024*1024), b''):h.update(part)
    return h.hexdigest()


def run(solver, case, output, seconds=240):
    solver=Path(solver).resolve();case=Path(case).resolve();output=Path(output)
    if not solver.is_file() or not case.is_file() or seconds<=0:
        raise ValueError('valid solver, case and positive timeout required')
    output.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    env.update(CFDX_V80_NATIVE_FACE_DEFECT='1',CFDX_SST_DIAGNOSTICS='1',
               CFDX_BFS_STAGE_TRACE='1')
    cmd=[str(solver),'--mesh',str(case),'--experimental-sst',
         '--iterations','1','--output-dir',str(output)]
    stderr=output/'native_first_sst.stderr';stdout=output/'native_first_sst.stdout'
    t=time.monotonic();reached=False;process_exit=None
    with stderr.open('w') as err,stdout.open('w') as out:
        p=subprocess.Popen(cmd,stdout=out,stderr=err,env=env)
        try:
            while time.monotonic()-t<seconds:
                if stderr.stat().st_size:
                    text=stderr.read_text(errors='replace')
                    if any(l.startswith('BFS_PHASE_V10') and
                           'phase=SST_TRANSPORT_END' in l for l in text.splitlines()):
                        reached=True;break
                if p.poll() is not None:break
                time.sleep(.25)
        finally:
            if p.poll() is None:
                p.kill() # stop before subsequent nonlinear retries
            process_exit=p.wait()
    status={
        'version':'V82-first-private-sst', 'solver_sha256':sha256(solver),
        'case_sha256':sha256(case),'first_sst_end_observed':reached,
        'process_exit_after_diagnostic_stop':process_exit,
        'elapsed_s':round(time.monotonic()-t,2),
        'outer_convergence_qualified':False,
        'reason':'bounded first private transaction; NOT a physical solution',
        'stderr':str(stderr),'stdout':str(stdout)}
    (output/'diagnostic_manifest_v82.json').write_text(json.dumps(status,indent=2))
    print(json.dumps(status,indent=2))
    if not reached:raise RuntimeError('first native SST transaction incomplete')
    return status


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--solver',required=True)
    p.add_argument('--case',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--timeout',type=int,default=240)
    args=p.parse_args()
    run(args.solver,args.case,args.out,args.timeout)

if __name__=='__main__':main()