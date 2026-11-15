#!/usr/bin/env python3
"""Enable the V88 local momentum step trust guard on a V87 BFS20540 HDF5.

Keep *every* mesh dataset, physical input, boundary patch, and final
convergence threshold byte-equivalent; change only the opt-in raw flag.
Never overwrite the user's original input.
"""
from pathlib import Path
import argparse, hashlib, json, shutil
import h5py
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(preimage, destination):
    preimage = Path(preimage).resolve(strict=True)
    destination = Path(destination).absolute()
    if destination == preimage or destination.exists():
        raise ValueError('Output must be new; never overwrite original')
    with h5py.File(preimage, 'r') as file:
        if int(file.attrs.get('n_cells', 0)) != 20540:
            raise ValueError('Requires BFS20540 mesh')
        case = json.loads(file.attrs['case_setup_json'])
        raw = case['numerics']['raw_settings']
        for requirement in ('v87_bfs_early_velocity_screen',
                            'v86_bfs_physical_trust_region',
                            'v85_bounded_nonmonotone_outer_progress'):
            if raw.get(requirement) != 'enabled':
                raise ValueError(f'V88 requires {requirement}=enabled')
        if 'v88_bfs_local_step_trust' in raw:
            raise ValueError('V88 option already exists')
        if case['numerics']['max_iterations'] != 100:
            raise ValueError('Requires same 100 outer iterations')
    shutil.copy2(preimage, destination)
    try:
        with h5py.File(destination, 'r+') as file:
            case = json.loads(file.attrs['case_setup_json'])
            case['numerics']['raw_settings']['v88_bfs_local_step_trust'] = 'enabled'
            file.attrs.modify('case_setup_json', json.dumps(case, separators=(',', ':')))
        with h5py.File(preimage, 'r') as before, h5py.File(destination, 'r') as after:
            if set(before.keys()) != set(after.keys()):
                raise RuntimeError('HDF5 groups changed')
            for name in before.attrs.keys():
                if name == 'case_setup_json':
                    old = json.loads(before.attrs[name]); new = json.loads(after.attrs[name])
                    del new['numerics']['raw_settings']['v88_bfs_local_step_trust']
                    if old != new:
                        raise RuntimeError('Case setup changed beyond V88 flag')
                elif name not in after.attrs or not np.array_equal(before.attrs[name], after.attrs[name]):
                    raise RuntimeError('Unexpected attribute modified: '+name)
            def check_dataset(name, obj):
                if isinstance(obj, h5py.Dataset):
                    if name not in after or not np.array_equal(obj[...], after[name][...]):
                        raise RuntimeError('Mesh/field dataset changed: '+name)
            before.visititems(check_dataset)
    except BaseException:
        destination.unlink(missing_ok=True)
        raise
    return {'preimage_sha256': sha(preimage), 'output_sha256': sha(destination),
            'cells': 20540, 'outer_target': 100,
            'only_change': 'v88_bfs_local_step_trust=enabled',
            'physical_validation': False}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('preimage', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.preimage, args.output), indent=2))