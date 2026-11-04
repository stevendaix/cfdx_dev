"""V22 diagnostics must not disguise non-convergence or extrapolation error."""
import numpy as np
import pytest

from scripts.diagnose_th03_v22 import (
    quadratic_wall_derivative, periodicity, momentum_hotspots, compare_wall_reconstructions,
)
from scripts.qualify_th03_v21 import H


def test_quadratic_wall_derivative_exact_for_quadratic_and_affine_fields():
    y1,y2=1.2e-5,4.8e-5
    for a,b in ((3.0e4,-1.4e8),(-5.0e3,0.0)):
        u1=a*y1+b*y1*y1
        u2=a*y2+b*y2*y2
        assert quadratic_wall_derivative(y1,u1,y2,u2)==pytest.approx(a,rel=1e-12)


def test_quadratic_rejects_illconditioned_and_nonfinite_stencils():
    for data in ((0.,1.,0.001,1.),(1.,1.,1.2,2.),(1.,1.,9.,2.),(1.,float('nan'),3.,2.)):
        with pytest.raises(ValueError): quadratic_wall_derivative(*data)


def test_limit_cycle_detects_period_six_without_claiming_physical_convergence():
    pattern=[.0058,.0290,.0125,.00975,.0289,.01251]
    history=[{'momentum_equation_residual_rhs_relative':pattern[i%6]} for i in range(120)]
    p=periodicity(history)
    assert p['best_lag']==6
    assert p['status']=='repeating_cycle_suspected'
    assert p['relative_rms_periodicity_error']<1e-12


def test_momentum_hotspot_coordinates_and_bad_index():
    centres=np.array([[-3*H, 2*H, 0], [2*H, H, 0]],float)
    bounds=np.array([[-3.1*H,-2.9*H,1.9*H,2.1*H],
                      [1.9*H,2.1*H,.9*H,1.1*H]])
    h=[{'momentum_residual_cell':0},{'momentum_residual_cell':0},{'momentum_residual_cell':1}]
    result=momentum_hotspots(h,centres,bounds)
    assert result['upstream_fraction']==pytest.approx(2/3)
    assert result['top'][0]['x_over_h']==pytest.approx(-3)
    with pytest.raises(ValueError): momentum_hotspots([{'momentum_residual_cell':2}],centres,bounds)


def test_v22_transport_scheme_opt_in_and_rejected_invalid_selection():
    from scripts.generate_th03_openfoam_blockmesh import make_case_setup
    def scheme(name):
        setup=make_case_setup(2,0.55,0.3,convection_scheme=name)
        return next(x.configuration_key for x in setup.numerics.selection.entries
                    if x.family=='convection')
    assert scheme('upwind')=='numerics.convection.upwind'
    assert scheme('second_order_upwind')=='numerics.convection.second_order_upwind'
    assert scheme('tvd_minmod')=='numerics.convection.tvd.minmod'
    with pytest.raises(ValueError,match='unsupported TH-03 convection'):
        make_case_setup(2,0.55,0.3,convection_scheme='central_unbounded')


def test_v22_tvd_mode_is_preserved_in_real_hdf5_case(tmp_path):
    import json
    import h5py
    from scripts.generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
    path=tmp_path/'th03_v22.cfdx.h5'
    result=write_case(path,make_structured_mesh('smoke'),mesh_source='V22 numerical selection test',
                      max_iterations=2,convection_scheme='tvd_minmod')
    assert result['momentum_convection_scheme']=='numerics.convection.tvd.minmod'
    with h5py.File(path,'r') as handle:
        raw=json.loads(handle.attrs['case_setup_json'])['numerics']['selection']['entries']
    assert any(x['family']=='convection' and
               x['configuration_key']=='numerics.convection.tvd.minmod' for x in raw)
