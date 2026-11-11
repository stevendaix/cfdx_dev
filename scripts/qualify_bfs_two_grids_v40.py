"""Fail-closed synthesis of BFS 712/20,540 six-Schur campaigns.

Three distinct evidence levels: attempted execution; algebraic/physical
convergence; independent Driver-Seegmiller agreement. No timeout or partial
solve is promoted to a physical pass. The supplied input is the *campaign
report*, not a saved marketing verdict. FNV64 is used only for within-run
operator identity and is not cryptographically collision resistant.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import statistics
from bfs_schur_fingerprint_v40 import compare_six_fingerprints
from generate_th03_bfs_case import SCHUR_MODELS

SCHEMA = 'CFDX-V40-BFS-2GRID-SCHUR-QUALIFICATION/FAIL-CLOSED'
LINEAR = re.compile(r'COUPLED_BLOCK_TRUE_RESIDUAL_V16[^\n]*?full_relative=([^\s]+)')


def read_pcd_csv(path: Path, n_cells: int) -> dict:
    with Path(path).open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    expected = {'cell', 'fp_scaled_shape_defect', 'kp_scaled_shape_defect',
                'gauge_like', 'fp_best_multiplier'}
    if len(rows) != n_cells or not rows or not expected <= set(rows[0]):
        raise ValueError('V40 PCD diagnostic incomplete or missing rows')
    values=[]; ids=set();worst=None
    for row in rows:
        i=int(row['cell']); v=float(row['fp_scaled_shape_defect'])
        if i<0 or i>=n_cells or i in ids or not math.isfinite(v) or not 0<=v<=1:
            raise ValueError('V40 PCD audit invalid cell/shape defect')
        ids.add(i);values.append(v)
        if worst is None or v>worst[1]:worst=(i,v)
    values.sort()
    return {'cells':len(rows),'median_shape_defect':statistics.median(values),
            'p90_shape_defect':values[int(.9*(len(values)-1))],
            'worst_cell':worst[0],'worst_shape_defect':worst[1],
            'meaning':'scale-invariant Sp_diag-vs-Fp stencil fit, NOT exact Schur/physical error'}


def summarize_one(data:dict,n_cells:int, pc_csv:Path|None=None) -> dict:
    if (not data.get('identical_physical_inputs') or
        int(next(iter(data.get('cases', {}).values()), {}).get('n_cells', -1)) != n_cells):
        raise ValueError('V40 malformed or wrong-size BFS campaign')
    results=data.get('results',[])
    if len(results)!=6 or {x.get('schur') for x in results}!=set(SCHUR_MODELS):
        raise ValueError('V40 missing or duplicate Schur variant')
    if len({x.get('physical_sha256') for x in data['cases'].values()})!=1:
        raise ValueError('V40 physical input mismatch')
    same=compare_six_fingerprints(results,SCHUR_MODELS)
    table=[]
    for r in results:
        if r['schur'] not in data['cases']:
            raise ValueError('V40 result does not match input Schur')
        completed = r.get('execution') in ('EXIT_ZERO','FAILED')
        physical = bool(r.get('nonlinear_convergence_proven'))
        if physical and (not completed or r.get('exit_code') != 0):
            raise ValueError('V40 invalid converged claim on incomplete/failed run')
        table.append({'schur':r['schur'], 'execution':r.get('execution'),
                      'first_blocker':r.get('first_blocker',{}).get('first_blocker'),
                      'elapsed_seconds':r.get('wall_seconds'),
                      'system_fingerprint_present':bool(r.get('operator_fingerprint_v40')),
                      'outer_physical_converged':physical,
                      'physical_validation':'NOT_QUALIFIED'})
    result={'cells':n_cells,'n_faces':next(iter(data['cases'].values()))['n_faces'],
            'identical_physical_hdf5':True,
            'all_six_same_frozen_A_b_x0_FNV64':same['all_six_same_A_b_x0_fingerprint'],
            'operator_proof_level':'FNV64 numerical equality, NON CRYPTOGRAPHIC',
            'full_matrix_capture_equal':data.get('all_six_identical_original_operator') is True,
            'schur_results':table,'all_six_outer_converged':all(x['outer_physical_converged'] for x in table),
            'driver_seegmiller_validated':False,
            'source_kind':data.get('source_kind')}
    if pc_csv is not None:
        result['pcd_schur_row_fit']=read_pcd_csv(pc_csv,n_cells)
    return result


def compile_report(first:dict,second:dict, pc_first=None, pc_second=None) -> dict:
    small=summarize_one(first,712,pc_first)
    large=summarize_one(second,20540,pc_second)
    report={'schema':SCHEMA,'grid_712':small,'grid_20540':large,
            'same_physical_input_across_different_grids':False,
            'grid_convergence_index_qualified':False,
            'global_U_p_k_omega_jacobian_qualified':False,
            'nasa_validation_qualified':False,
            'all_twelve_outer_converged':small['all_six_outer_converged'] and large['all_six_outer_converged'],
            'limits':['reconstructed graded six-block mesh, not native OpenFOAM polyMesh',
                      'no independent NASA fields, reattachment or GCI',
                      'incomplete runs cannot certify either linear or nonlinear convergence']}
    json.dumps(report,allow_nan=False)
    return report


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign-712',type=Path,required=True)
    p.add_argument('--campaign-20540',type=Path,required=True)
    p.add_argument('--csv-712',type=Path)
    p.add_argument('--csv-20540',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(argv)
    try:
        d1=json.loads(a.campaign_712.read_text())
        d2=json.loads(a.campaign_20540.read_text())
        report=compile_report(d1,d2,a.csv_712,a.csv_20540)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        tmp=a.output.with_name(a.output.name+'.partial')
        tmp.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
        os.replace(tmp,a.output)
        ok=report['all_twelve_outer_converged'] and report['nasa_validation_qualified']
        print('V40_BFS_SIX_SCHUR_2GRIDS attempted=12 frozen_same_712={} frozen_same_20540={} '
              'all_outer_converged={} nasa_validated={} qualification={}'.format(
               report['grid_712']['all_six_same_frozen_A_b_x0_FNV64'],
               report['grid_20540']['all_six_same_frozen_A_b_x0_FNV64'],
               report['all_twelve_outer_converged'], report['nasa_validation_qualified'],
               'PASS' if ok else 'FAIL_OPEN_GATES'))
        return 0 if ok else 2
    except (ValueError,KeyError,TypeError,OSError,json.JSONDecodeError) as exc:
        p.exit(3, f'V40 BFS qualification cannot be evaluated: {exc}\n')

if __name__=='__main__':raise SystemExit(main())
