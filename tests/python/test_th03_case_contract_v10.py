"""V10: reject false TH-03 input equivalence and reconstruct inlet omega."""
from __future__ import annotations

import json
import math
from pathlib import Path

import h5py
import pytest

from scripts import generate_th03_openfoam_blockmesh as g
from scripts.audit_th03_case_contract_v10 import audit


def _new_case(path: Path, *, omega: float = g.OMEGA_IN,
              omega_initial: float = g.OMEGA_INITIAL) -> Path:
    g.write_case(path, g.make_structured_mesh('smoke'),
                 mesh_source='independent six-block TH-03 fixture',
                 inlet_omega=omega, initial_omega=omega_initial,
                 max_iterations=2)
    return path


def test_documented_turbulent_viscosity_ratio_is_respected():
    assert g.K_IN / (g.NU * g.OMEGA_IN) == pytest.approx(.009, rel=1e-14)
    assert g.OMEGA_IN == pytest.approx(7763.532763532765)
    assert g.K_IN / (g.NU * g.OMEGA_LEGACY_V9) == pytest.approx(.00038448558)


def test_initial_omega_and_inlet_omega_are_independent():
    c = g.make_case_setup(2, .25, .1, initial_omega=g.OMEGA_LEGACY_V9)
    bc = next(x for x in c.boundary_conditions if x.patch_name == 'inlet')
    inlet_omega = next(x.value for x in bc.sst_scalar_conditions if x.field == 'omega')
    assert inlet_omega == pytest.approx(g.OMEGA_IN)
    assert c.initial_condition.scalar_fields['omega'] == pytest.approx(g.OMEGA_LEGACY_V9)


@pytest.mark.parametrize('bad', [float('nan'), float('inf'), -1., 0.])
def test_invalid_inlet_or_initial_omega_rejected(bad: float):
    with pytest.raises(ValueError, match='omega'):
        g.make_case_setup(2, .25, .1, inlet_omega=bad)
    with pytest.raises(ValueError, match='omega'):
        g.make_case_setup(2, .25, .1, initial_omega=bad)


def test_preflight_separates_reconstructed_openfoam_and_experiment(tmp_path: Path):
    case = _new_case(tmp_path / 'nominal.h5')
    result = audit(case)
    assert result['cells'] == 712
    assert result['mesh_status'] == 'SURROGATE_GEOMETRY_ONLY'
    assert result['physical_qualification'] == 'NOT_QUALIFIED'
    assert result['experimental_comparability'] == 'NOT_COMPARABLE'
    assert result['checks']['inlet_nut_over_nu_matches_documented_openfoam_0p009']
    assert result['checks']['inlet_turbulence_intensity_matches_nominal_openfoam_0p061percent']
    assert result['physics']['inlet_turbulence_intensity'] == pytest.approx(0.0006095, rel=0.01)
    assert not result['checks']['experimental_Re_h_37500']
    assert result['checks']['straight_recovery_floor_at_outlet']
    assert result['dimensions']['outlet_lower_wall_y_over_h'] == pytest.approx(0.0)
    assert result['physics']['reynolds_h'] == pytest.approx(35983.33333333)
    assert result['physics']['estimated_first_cell_yplus_from_upstream_experimental_cf'] > 0
    assert result['physics']['wall_resolved_sst']
    assert not result['checks']['second_order_convection_selected']


def test_preflight_diagnoses_legacy_inlet_mismatch(tmp_path: Path):
    case = _new_case(tmp_path / 'legacy.h5', omega=g.OMEGA_LEGACY_V9)
    report = audit(case)
    assert not report['checks']['inlet_nut_over_nu_matches_documented_openfoam_0p009']
    assert report['physics']['inlet_nut_over_nu'] == pytest.approx(.00038448558)
    assert report['physical_qualification'] == 'NOT_QUALIFIED'


def test_preflight_rejects_missing_explicit_k_omega(tmp_path: Path):
    case = _new_case(tmp_path / 'invalid.h5')
    with h5py.File(case, 'r+') as h:
        setup = json.loads(h.attrs['case_setup_json'])
        inlet = next(x for x in setup['boundary_conditions'] if x['patch_name'] == 'inlet')
        inlet['sst_scalar_conditions'] = [x for x in inlet['sst_scalar_conditions'] if x['field'] != 'omega']
        h.attrs['case_setup_json'] = json.dumps(setup)
    with pytest.raises(ValueError, match='omega'):
        audit(case)


def test_geometry_preflight_rejects_face_inversion_with_unchanged_bounds(tmp_path: Path):
    import numpy as np
    case = _new_case(tmp_path / 'inverted.h5')
    with h5py.File(case, 'r+') as h:
        neighbours = h['neighbour'][:]
        face = int(np.flatnonzero(neighbours >= 0)[12])
        offsets = h['face_offsets'][:]
        first, last = map(int, offsets[face:face+2])
        faces = h['face_vertices'][:]
        faces[first:last] = faces[first:last][::-1]
        h['face_vertices'][:] = faces
    report = audit(case)
    assert report['mesh_status'] == 'FAIL_GEOMETRY'
    assert report['mesh_quality']['flipped_internal_faces'] >= 1
    assert report['physical_qualification'] == 'NOT_QUALIFIED'


def test_re_is_not_silently_changed_to_match_experiment(tmp_path: Path):
    case = _new_case(tmp_path / 'original.h5')
    original = audit(case)
    assert not original['checks']['experimental_Re_h_37500']
    with h5py.File(case, 'r+') as h:
        setup = json.loads(h.attrs['case_setup_json'])
        setup['materials'][0]['dynamic_viscosity'] = g.RHO * g.VELOCITY * g.H / 37500
        h.attrs['case_setup_json'] = json.dumps(setup)
    changed = audit(case)
    assert changed['checks']['experimental_Re_h_37500']
    assert changed['experimental_comparability'] == 'NOT_COMPARABLE'
    assert not changed['checks']['measured_x_over_h_minus4_profile_declared']
