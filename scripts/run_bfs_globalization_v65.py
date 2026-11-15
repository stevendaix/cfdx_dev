"""Prepare (or explicitly execute) V65 BFS152 three-way comparison.

Three executions share one HDF5 geometry/physics/BC/SST setup. The only
experimental changes are V13/V14/V47 and the V65 state-only merit. This is a
DIAGNOSTIC campaign, never an assertion of NASA/OpenFOAM qualification.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

try:
    import h5py
except ImportError as ex:
    raise SystemExit('V65 requires h5py for a reproducible HDF5 case') from ex

from audit_bfs_globalization_v65 import audit,physical_setup_equal

FLAGS={
    'v13_joint_globalization':'true',
    'v14_monotone_global_merit':'true',
    'coupled_composite_backtracking_v47_experimental':'true',
    'coupled_composite_backtracking_retries_v47':'10',
}


def prepare_case(src:Path,dest:Path,steps:int,enable_v47:bool):
    if steps<1 or steps>200:raise ValueError('bounded BFS152 campaign requires 1..200 steps')
    shutil.copy2(src,dest)
    with h5py.File(dest,'r+') as f:
        if int(f.attrs['n_cells'])!=152 or int(f.attrs['n_faces'])!=634:
            raise ValueError('V65 currently only qualified for BFS152 topology')
        settings=json.loads(f.attrs['case_setup_json'])
        if settings['numerics'].get('pressure_velocity') not in (None,'coupled'):
            # Most CFDX cases store this selection in raw_settings; leave
            # resolution of other spellings to the native production parser.
            pass
        settings['numerics']['max_iterations']=steps
        raw=settings['numerics']['raw_settings']
        for flag in FLAGS:raw.pop(flag,None)
        if enable_v47:raw.update(FLAGS)
        f.attrs['case_setup_json']=json.dumps(settings,separators=(',',':'))


def run_one(solver:Path,case:Path,dest:Path,deadline:float,v65:bool):
    dest.mkdir(exist_ok=False)
    env={k:v for k,v in os.environ.items()
         if not (k.startswith('CFDX_FPHYS_') or k.startswith('CFDX_POST_SST_') or
                 k.startswith('CFDX_BFS_') or k.startswith('CFDX_COUPLED_SELFP_'))}
    if v65:
        env['CFDX_FPHYS_CANDIDATE_V51']='1'
        env['CFDX_FPHYS_V65_SAME_OPERATOR']='1'
    env['CFDX_BFS_STAGE_TRACE']='1'
    args=[str(solver),'--mesh',str(case),'--experimental-sst',
          '--output-dir',str(dest/'result')]
    with (dest/'solver.log').open('w') as output:
        try:
            result=subprocess.run(args,env=env,stdout=output,stderr=subprocess.STDOUT,
                                  timeout=deadline,check=False)
            code=result.returncode
        except subprocess.TimeoutExpired:
            code=124
    return code


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base-case',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--solver',type=Path)
    ap.add_argument('--iterations',type=int,default=40)
    ap.add_argument('--wall-timeout',type=float,default=900.)
    ap.add_argument('--execute',action='store_true',help='Opt-in three long native solver runs')
    ap.add_argument('--allow-partial',action='store_true',help='Only for an explicitly incomplete diagnostic')
    args=ap.parse_args()
    if args.out.exists():ap.error('refusing to overwrite evidence directory')
    original=args.base_case.resolve(strict=True)
    args.out.mkdir(parents=True)
    baseline=args.out/'bfs152_baseline.cfdx.h5'
    v47=args.out/'bfs152_v47.cfdx.h5'
    prepare_case(original,baseline,args.iterations,False)
    prepare_case(original,v47,args.iterations,True)
    physical_setup_equal(baseline,v47)
    print('V65_CASES_PREPARED identical_physics=1 cells=152 steps=',args.iterations)
    if not args.execute:return
    if not args.solver:ap.error('--execute requires --solver')
    solver=args.solver.resolve(strict=True)
    results={}
    for name,case,v65 in [('coupled',baseline,False),('v47_old',v47,False),
                           ('v47_v65',v47,True)]:
        results[name]=run_one(solver,case,args.out/name,args.wall_timeout,v65)
    log=lambda name:(args.out/name/'solver.log').read_text(errors='replace')
    result=audit(log('v47_v65'),log('coupled'),log('v47_old'))
    result['cases']={'coupled':str(baseline),'v47':str(v47)}
    result['exit_codes']=results
    result['solver_sha256']=hashlib.sha256(solver.read_bytes()).hexdigest()
    (args.out/'audit.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    if results['v47_v65']==124 or not result['complete_40_iteration_run']:
        if not args.allow_partial:
            raise RuntimeError('V65 long campaign incomplete: see audit.json, never a PASS')
    print('V65_COMPARISON_RECORDED steps=',result['iterations_committed'],
          ' qualified=0')

if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError) as ex:
        print('V65_RUN_REJECTED:',str(ex),file=sys.stderr)
        raise SystemExit(2)