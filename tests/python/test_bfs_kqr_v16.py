"""Physical BC identity: OpenFOAM kqR wall k is zero-gradient, never k=0."""
import json
from pathlib import Path
import h5py
from scripts.generate_th03_openfoam_blockmesh import make_case_setup, make_structured_mesh, write_case


def test_opt_in_kqr_wall_in_setup():
    old = make_case_setup(2, .25, .1)
    new = make_case_setup(2, .25, .1, openfoam_kqr_wall_v16=True)
    for setup, expected in ((old, 'wall_k_zero'),(new,'wall_kq_r')):
        for patch in setup.boundary_conditions:
            if patch.type.value != 'wall':
                continue
            sc = {a.field: a.type.value for a in patch.sst_scalar_conditions}
            assert sc == {'k': expected, 'omega': 'wall_omega_sst'}


def test_kqr_survives_hdf5(tmp_path: Path):
    path = tmp_path/'bfs.cfdx.h5'
    write_case(path, make_structured_mesh('smoke'),mesh_source='geometric-surrogate',
               max_iterations=1,openfoam_kqr_wall_v16=True)
    with h5py.File(path) as f:
        setup=json.loads(f.attrs['case_setup_json'])
        for patch in setup['boundary_conditions']:
            if patch['type'] == 'wall':
                assert {s['field']:s['type'] for s in patch['sst_scalar_conditions']} == {
                    'k':'wall_kq_r','omega':'wall_omega_sst'}
