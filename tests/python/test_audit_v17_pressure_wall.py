"""Fail-closed evidence parser: flux topology vs physical qualification."""
from pathlib import Path
import pytest
from scripts.audit_v17_pressure_wall import audit_log


def sample_log() -> str:
    return ('V17_FLUX_AUDIT iteration=1 cell_l1=0.2 cell_linf=0.2 '
            'cell_worst=2 worst_div=0.2 boundary_signed=0.2 boundary_abs=2.2 '
            'cell_signed=0.2 topology_error=0 rhie_chow_identity_linf=0 qualified=0\n'
            'V17_FLUX_PATCH iteration=1 name=inlet faces=1 signed=-1 absolute=1 qualified=0\n'
            'V17_FLUX_PATCH iteration=1 name=outlet faces=1 signed=1.2 absolute=1.2 qualified=0\n')


def write(tmp_path: Path, text: str) -> Path:
    p=tmp_path/'log.txt'
    p.write_text(text)
    return p


def test_v17_audit_closes_by_patch_and_refuses_to_qualify(tmp_path):
    r=audit_log(write(tmp_path,sample_log()))
    assert r['final_iteration'] == 1
    assert r['status'] == 'DIAGNOSTIC_ONLY_NOT_QUALIFIED'
    assert r['wall']['physical_yplus_qualified'] is False


def test_v17_audit_rejects_wrong_boundary_signed(tmp_path):
    p=write(tmp_path,sample_log().replace('signed=1.2','signed=1.1'))
    with pytest.raises(ValueError,match='patch signed'):
        audit_log(p)


def test_v17_audit_rejects_missing_summary(tmp_path):
    with pytest.raises(ValueError, match='no V17'):
        audit_log(write(tmp_path,'Converged YES Iterations 1\n'))


def test_v17_audit_rejects_nonfinite_face(tmp_path):
    p=write(tmp_path,sample_log().replace('cell_linf=0.2','cell_linf=nan'))
    with pytest.raises(ValueError, match='nonfinite'):
        audit_log(p)


def test_v17_audit_checks_pseudo_mass(tmp_path):
    extra = ('V17_PSEUDO_MOMENTUM iteration=1 cfl=0.3 mass_min=0.1 mass_max=2 '
             'physical_momentum_rhs_relative=0.15 continuity=0.01 qualified=0\n')
    result=audit_log(write(tmp_path,sample_log()+extra))
    assert result['pseudo'][0]['cfl']==0.3


def test_v17_audit_refuses_pseudo_qualifier(tmp_path):
    extra = ('V17_PSEUDO_MOMENTUM iteration=1 cfl=0.3 mass_min=0.1 mass_max=2 '
             'physical_momentum_rhs_relative=0.15 continuity=0.01 qualified=1\n')
    with pytest.raises(ValueError,match='cannot be marked qualified'):
        audit_log(write(tmp_path,sample_log()+extra))


def test_v17_pressure_correction_audit_is_independent_of_linear_residual(tmp_path):
    extra = ('V17_PRESSURE_CORRECTION iteration=1 corrector=0 alpha_p=0.1 '
             'before_l1=2 after_l1=1.5 before_linf=1 after_linf=0.8 '
             'reduction_ratio=0.75 pressure_linear_relative=1e-12 qualified=0\n')
    r=audit_log(write(tmp_path,sample_log()+extra))
    assert r['pressure_corrections'][0]['reduction_ratio']==0.75
    assert r['status'] == 'DIAGNOSTIC_ONLY_NOT_QUALIFIED'


def test_v17_pressure_correction_rejects_mismatched_ratio(tmp_path):
    extra = ('V17_PRESSURE_CORRECTION iteration=1 corrector=0 alpha_p=0.1 '
             'before_l1=2 after_l1=1.5 before_linf=1 after_linf=0.8 '
             'reduction_ratio=0.5 pressure_linear_relative=1e-12 qualified=0\n')
    with pytest.raises(ValueError,match='ratio mismatch'):
        audit_log(write(tmp_path,sample_log()+extra))


def test_v17_worst_cell_face_flux_closes(tmp_path):
    extra = ('V17_FLUX_FACE iteration=1 cell=2 face=4 patch=internal outward=0.3 '
             'owner_flux=-0.3 qualified=0\n'
             'V17_FLUX_FACE iteration=1 cell=2 face=5 patch=outlet outward=-0.1 '
             'owner_flux=-0.1 qualified=0\n')
    r=audit_log(write(tmp_path,sample_log()+extra))
    assert len(r['final_flux']['worst_cell_faces'])==2


def test_v17_worst_cell_face_flux_rejects_nonclosure(tmp_path):
    extra = ('V17_FLUX_FACE iteration=1 cell=2 face=4 patch=internal outward=0.3 '
             'owner_flux=-0.3 qualified=0\n'
             'V17_FLUX_FACE iteration=1 cell=2 face=5 patch=outlet outward=0.1 '
             'owner_flux=0.1 qualified=0\n')
    with pytest.raises(ValueError,match='worst cell face fluxes'):
        audit_log(write(tmp_path,sample_log()+extra))