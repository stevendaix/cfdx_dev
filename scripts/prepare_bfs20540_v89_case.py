#!/usr/bin/env python3
"""Create opt-in V89 BFS20 540 case; verify exact HDF5 mesh and physics.

V89 adds ONE raw flag to an already prepared V88 case. No input overwritten.
"""
import argparse, hashlib, json, shutil
from pathlib import Path
import h5py
import numpy as np

def prepare(preimage, destination):
    preimage=Path(preimage).resolve(strict=True)
    destination=Path(destination).absolute()
    if destination==preimage or destination.exists():
        raise ValueError('Source may not be overwritten')
    with h5py.File(preimage) as f:
        if int(f.attrs['n_cells'])!=20540:raise ValueError('Only BFS20 540')
        case=json.loads(f.attrs['case_setup_json'])
        raw=case['numerics']['raw_settings']
        if raw.get('v88_bfs_local_step_trust')!='enabled':
            raise ValueError('V89 requires the V88 exact preimage')
        if 'v89_bfs_pressure_transverse_trust' in raw:
            raise ValueError('V89 already set')
    shutil.copy2(preimage,destination)
    try:
        with h5py.File(destination,'r+') as f:
            case=json.loads(f.attrs['case_setup_json'])
            case['numerics']['raw_settings']['v89_bfs_pressure_transverse_trust']='enabled'
            f.attrs.modify('case_setup_json',json.dumps(case,separators=(',',':')))
        with h5py.File(preimage) as a,h5py.File(destination) as b:
            if set(a.keys())!=set(b.keys()) or set(a.attrs)!=set(b.attrs):
                raise ValueError('HDF5 topology modified')
            def check(name,obj):
                if isinstance(obj,h5py.Dataset):
                    if name not in b or not np.array_equal(obj[...],b[name][...]):
                        raise ValueError('Changed dataset: '+name)
            a.visititems(check)
            for name in a.attrs:
                if name=='case_setup_json':
                    old=json.loads(a.attrs[name]); new=json.loads(b.attrs[name]);
                    del new['numerics']['raw_settings']['v89_bfs_pressure_transverse_trust']
                    if old!=new: raise ValueError('Other numerical settings modified')
                elif not np.array_equal(a.attrs[name],b.attrs[name]):
                    raise ValueError('Changed attribute '+name)
    except Exception:
        destination.unlink(missing_ok=True)
        raise
    return {'n_cells':20540,'only_change':'v89_bfs_pressure_transverse_trust=enabled',
            'source_sha256':hashlib.sha256(preimage.read_bytes()).hexdigest(),
            'output_sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),
            'qualified':False}
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('preimage',type=Path)
    ap.add_argument('output',type=Path)
    a=ap.parse_args()
    print(json.dumps(prepare(a.preimage,a.output),indent=2))