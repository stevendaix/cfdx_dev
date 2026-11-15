"""Prepare non-destructive 100-outer-iteration V85 six-block BFS case.

Never changes geometry, material, initial conditions, turbulent wall model,
solver equations, or final physical convergence tolerances.
"""
from pathlib import Path
import argparse, hashlib, json, shutil
import h5py

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def prepare(source,destination,count=100):
    source,destination=Path(source),Path(destination)
    if source.resolve()==destination.resolve():raise ValueError('Refuse to overwrite reference case')
    if not (1<=count<=100):raise ValueError('Campaign cap must be 1..100')
    with h5py.File(source,'r') as f:
        if int(f.attrs['n_cells'])!=20540 or str(f.attrs.get('source_case_name',''))!='backwardFacingStep2D':
            raise ValueError('Not the six-block BFS20 540 reference geometry')
        setup=json.loads(f.attrs['case_setup_json'])
        raw=setup['numerics']['raw_settings']
        if setup['numerics']['coupled_solver']!='SIMPLE' or setup.get('transient',False):
            raise ValueError('V85 only for steady SIMPLE SST')
        for flag in ('v7_post_sst_momentum_gate','v8_global_nonlinear_residual',
                     'v13_joint_globalization','v14_monotone_global_merit',
                     'pressure_high_low_production'):
            if raw.get(flag)!='true':raise ValueError('Missing strict baseline gate '+flag)
        if float(raw.get('pressure_linear_tolerance','1'))>1e-11:
            raise ValueError('Pressure tolerance must remain at or below 1e-11')
        if int(raw.get('sst_transport_max_iterations','0'))<10:
            raise ValueError('SST budget must be at least 10')
    destination.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,destination)
    with h5py.File(destination,'r+') as f:
        setup=json.loads(f.attrs['case_setup_json']);raw=setup['numerics']['raw_settings']
        setup['numerics']['max_iterations']=count
        raw['sst_outer_provisional']='enabled'
        raw['v83_provisional_outer_merit_progress']='enabled'
        raw['v84_bounded_nonmonotone_outer_progress']='disabled'
        raw['v85_bounded_nonmonotone_outer_progress']='enabled'
        f.attrs.modify('case_setup_json',json.dumps(setup,separators=(',',':')))
    with h5py.File(source,'r') as a,h5py.File(destination,'r') as b:
        for h in ('mesh_hash','topology_hash'):
            if a.attrs.get(h)!=b.attrs.get(h):raise RuntimeError('Mesh hash changed '+h)
        # Physical datasets must remain byte-identical, including connectivity.
        def compare(name,obj):
            if isinstance(obj,h5py.Dataset):
                if name not in b:raise RuntimeError('Missing dataset '+name)
                import numpy as np
                if not np.array_equal(obj[...],b[name][...],equal_nan=True):
                    raise RuntimeError('Case dataset mutated: '+name)
        a.visititems(compare)
    return {'source':str(source),'output':str(destination),'n_cells':20540,
      'source_sha256':digest(source),'output_sha256':digest(destination),
      'max_outer_iterations':count,'vtu_stride':10,'final_qualified':False,
      'changes':['numerics.max_iterations','numerics.raw_settings.sst_outer_provisional',
                 'numerics.raw_settings.v83_provisional_outer_merit_progress',
                 'numerics.raw_settings.v84_bounded_nonmonotone_outer_progress',
                 'numerics.raw_settings.v85_bounded_nonmonotone_outer_progress']}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('destination',type=Path)
    parser.add_argument('--max-outer',type=int,default=100)
    a=parser.parse_args();print(json.dumps(prepare(a.source,a.destination,a.max_outer),indent=2))