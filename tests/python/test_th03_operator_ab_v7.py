"""The face-policy campaign preserves the strict pressure and physical gates."""
from pathlib import Path
import json
import h5py
import pytest

from scripts.run_th03_operator_ab_v7 import compare


def test_rejects_invalid_controls(tmp_path):
    for bad in (0.,-1.,float('nan')):
        with pytest.raises(ValueError):
            compare(tmp_path/'bad',tmp_path/'run', ('smoke',),timeout=bad)
    with pytest.raises(ValueError):
        compare(tmp_path/'bad',tmp_path/'run', ('smoke',),picard_damping=1.5)


def test_run_records_both_policies_but_never_fakes_physical_pass(tmp_path):
    fake = tmp_path / 'solver.sh'
    fake.write_text('#!/bin/sh\necho "solve_steady_incompressible: pressure-correction solve did not converge (status=1, iterations=13, residual=1e-12, relative=3e-11)"\nexit 2\n')
    fake.chmod(0o755)
    out = tmp_path / 'run'
    result = compare(fake,out,('smoke',),timeout=5)
    assert result['pressure_tolerance'] == 1e-11
    assert result['qualification'] == 'NOT_QUALIFIED'
    assert result['cf_cp_xr'] == 'UNAVAILABLE'
    assert len(result['entries']) == 2
    assert {entry['policy'] for entry in result['entries']} == {'arithmetic','projected'}
    assert all(entry['pressure_failure']['relative_residual']==3e-11
               for entry in result['entries'])
    for policy in ('arithmetic','projected'):
        with h5py.File(out/f'smoke_{policy}.cfdx.h5') as h5:
            s=json.loads(h5.attrs['case_setup_json'])
        assert s['numerics']['raw_settings']['pressure_face_interpolation']==policy
        assert s['numerics']['raw_settings']['pressure_linear_tolerance']=='1e-11'
    result_diag = compare(fake,tmp_path/'diag',('smoke',),diagnostic=True,
                          outer_provisional=True)
    assert result_diag['mode']=='DIAGNOSTIC_ONLY'
    assert result_diag['pressure_tolerance']==1e-9
    assert result_diag['qualification']=='NOT_QUALIFIED'


def test_provisional_mode_is_recorded_but_never_qualified(tmp_path):
    fake = tmp_path/'noop.sh'
    fake.write_text('#!/bin/sh\necho "SST_PROVISIONAL_OUTER iteration=1 qualified=0"\nexit 2\n')
    fake.chmod(0o755)
    result=compare(fake,tmp_path/'p',('smoke',),diagnostic=True,
                   picard_damping=.2,outer_provisional=True)
    assert result['qualification']=='NOT_QUALIFIED'
    assert result['sst_provisional_outer'] is True
    assert all(e['physical_qualification']=='NOT_QUALIFIED' for e in result['entries'])
    with h5py.File(tmp_path/'p/smoke_projected.cfdx.h5') as h5:
        raw=json.loads(h5.attrs['case_setup_json'])['numerics']['raw_settings']
        assert raw['sst_outer_provisional']=='enabled'
        assert raw['sst_picard_damping']=='0.2'
