"""A converged but physically wrong Cf/profile must never be certified."""
from pathlib import Path
import pytest
from scripts.qualify_bfs_nasa_fields_v16 import (
    read_nasa_cf,read_nasa_cp,read_nasa_profiles,compare_cf,compare_cp,compare_reattachment,
    compare_profiles,audit_data,ComparisonError)


def test_cf_uncertainty_finds_wrong_solution(tmp_path:Path):
    path=tmp_path/'cf.exp.dat'
    path.write_text('variables="x","cf","error"\n' +
        '\n'.join(f'{i} 0.001 0.0001' for i in range(12)))
    ref=read_nasa_cf(path)
    same=[{'x_h':x['x_h'],'cf':x['cf']} for x in ref]
    wrong=[{'x_h':x['x_h'],'cf':x['cf']+0.003} for x in ref]
    assert compare_cf(ref,same)['all_points_within_3_sigma']
    result=compare_cf(ref,wrong)
    assert not result['all_points_within_3_sigma']
    assert result['max_abs_experimental_sigma']>20


def test_no_experimental_data_may_not_pass(tmp_path:Path):
    with pytest.raises(ComparisonError,match='all required'):
        audit_data(tmp_path,{'cf':[],'profiles':{}})


def test_cp_offset_cannot_be_mistaken_for_original(tmp_path:Path):
    p=tmp_path/'cp.expnew.dat'
    p.write_text('zone,t="bottom wall"\n'+
        '\n'.join(f'{i} 0.1 0.2' for i in range(5))+
        '\nzone,t="top wall"\n'+
        '\n'.join(f'{i} 0.3 0.4' for i in range(5)))
    data=read_nasa_cp(p)
    assert data['bottom wall'][0]['cp_shifted']==0.1
    assert data['bottom wall'][0]['cp_original']==0.2


def test_missing_reynolds_stress_profile_is_rejected(tmp_path:Path):
    p=tmp_path/'profiles.exp.dat'
    p.write_text('variables="ID","y","u","v","uu","vv","uv","uuu","uvv","vuu","vvv"\n'+
        ''.join(f'zone,t="exp, x/H={station:g}"\n'+
                '\n'.join(f'{i} {i*.1} 0.5 0.2 1.0 1.0 -10.0 0 0 0 0' for i in range(1,6))+'\n'
                for station in (-4.,1.,4.,6.,10.)))
    ref=read_nasa_profiles(p)
    with pytest.raises(ComparisonError,match='missing modeled profile'):
        compare_profiles(ref,{})
    assert ref[1.][0]['uv_ref2']==-0.01


def test_cp_origin_mismatch_is_detected(tmp_path:Path):
    p=tmp_path/'cp.expnew.dat'
    p.write_text('zone,t="bottom wall"\n'+
        '\n'.join(f'{i} 0.1 0.2' for i in range(5))+
        '\nzone,t="top wall"\n'+
        '\n'.join(f'{i} 0.3 0.4' for i in range(5)))
    reference=read_nasa_cp(p)
    predicted={wall:[{'x_h':row['x_h'],'cp':row['cp_shifted']} for row in rows]
               for wall,rows in reference.items()}
    assert compare_cp(reference,predicted,'shifted')['walls']['bottom wall']['max_abs_error']==0.
    wrong=compare_cp(reference,predicted,'original')
    assert wrong['walls']['bottom wall']['max_abs_error']>=0.099
    with pytest.raises(ComparisonError,match='convention'):
        compare_cp(reference,predicted,'unspecified')


def test_reattachment_does_not_qualify_wrong_physics():
    good=compare_reattachment({'reattachment_x_h':6.27,
        'reattachment_method':'wall_shear_zero_crossing'})
    assert good['within_experimental_band'] and not good['qualified']
    wrong=compare_reattachment({'reattachment_x_h':5.3,
        'reattachment_method':'wall_shear_zero_crossing'})
    assert not wrong['within_experimental_band']
    with pytest.raises(ComparisonError,match='wall shear'):
        compare_reattachment({'reattachment_x_h':6.26,'reattachment_method':'eye_balling'})


def test_small_nonlinear_residual_cannot_hide_physically_wrong_solution(tmp_path:Path):
    """Synthetic oracle fixture, NOT actual NASA data; checks independent gates."""
    (tmp_path/'cf.exp.dat').write_text('\n'.join(
        f'{i} 0.001 0.0001' for i in range(12))+'\n')
    (tmp_path/'cp.expnew.dat').write_text(
        'zone,t="bottom wall"\n'+'\n'.join(f'{i} 0.1 0.2' for i in range(5))+
        '\nzone,t="top wall"\n'+'\n'.join(f'{i} 0.3 0.4' for i in range(5)))
    (tmp_path/'profiles.exp.dat').write_text(
        ''.join(f'zone,t="exp, x/H={station:g}"\n'+
                '\n'.join(f'{i} {i*.1} 0.5 0.2 1.0 1.0 -10.0 0 0 0 0'
                          for i in range(1,6))+'\n'
                for station in (-4.,1.,4.,6.,10.)))
    predicted={
        'nonlinear_relative_residual':1e-14,
        'cf_method':'independent_second_order_wall_shear',
        'cf':[{'x_h':i,'cf':0.006} for i in range(12)],
        'cp_convention':'original',
        'cp':{
            'bottom wall':[{'x_h':i,'cp':0.5} for i in range(5)],
            'top wall':[{'x_h':i,'cp':0.5} for i in range(5)]},
        'reattachment_x_h':5.3,
        'reattachment_method':'wall_shear_zero_crossing',
        'profiles':{str(int(station)):[
            {'y_h':i*.1,'u_ref':0.,'v_ref':0.,'uv_ref2':0.}
            for i in range(1,6)] for station in (-4.,1.,4.,6.,10.)}}
    evidence=audit_data(tmp_path,predicted)
    assert evidence['physical_validation']=='NOT_QUALIFIED'
    assert not evidence['reattachment_diagnostic']['within_experimental_band']
    assert not evidence['cf_diagnostic']['all_points_within_3_sigma']
    assert evidence['cp_diagnostic']['walls']['bottom wall']['max_abs_error']>=0.3-1e-12
    assert evidence['profile_diagnostic']['profiles']['1']['uv_ref2']['max_abs_error']>0
