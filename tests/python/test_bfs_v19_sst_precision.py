"""V19 opt-in SST precision fix and OpenFOAM corner boundary contract."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import generate_th03_openfoam_blockmesh as bfs


def test_v19_defaults_disabled_and_nasa_not_qualified(tmp_path: Path) -> None:
    mesh = bfs.make_structured_mesh('smoke')
    path = tmp_path / 'default.cfdx.h5'
    bfs.write_case(path, mesh, mesh_source='V19-smoke')
    with h5py.File(path) as h:
        data = json.loads(h.attrs['case_setup_json'])
    settings = data['numerics']['raw_settings']
    assert settings['v19_sst_equilibrated_defect_correction_experimental'] == 'disabled'
    assert settings['v19_omega_corner_equal_faces_experimental'] == 'disabled'
    assert data['source_metadata']['physical_qualification'] == 'NOT_QUALIFIED'


def test_v19_omega_corner_requires_v18_owners(tmp_path: Path) -> None:
    mesh = bfs.make_structured_mesh('smoke')
    with pytest.raises(ValueError, match='owner|V18|v18'):
        bfs.write_case(tmp_path/'invalid.cfdx.h5', mesh,
                       mesh_source='V19-smoke',v19_omega_corner_equal_faces=True)


def test_v19_two_modes_are_persisted_without_changing_sst_tolerance(tmp_path: Path) -> None:
    mesh = bfs.make_structured_mesh('smoke')
    path = tmp_path/'v19.cfdx.h5'
    bfs.write_case(path, mesh, mesh_source='V19-smoke',
                   sst_molecular_wall_diffusion=True,
                   openfoam_kqr_wall_v16=True,
                   v17_active_wall_operator=True,
                   v18_omega_owner_cells=True,
                   v19_sst_equilibrated_defect=True,
                   v19_omega_corner_equal_faces=True,
                   sst_inner_iterations=800)
    with h5py.File(path) as h:
        data = json.loads(h.attrs['case_setup_json'])
        assert int(h.attrs['n_cells']) == 712
    settings = data['numerics']['raw_settings']
    assert settings['v19_sst_equilibrated_defect_correction_experimental'] == 'enabled'
    assert settings['v19_omega_corner_equal_faces_experimental'] == 'enabled'
    assert settings['sst_transport_max_iterations'] == '800'
    assert data['source_metadata']['physical_qualification'] == 'NOT_QUALIFIED'
