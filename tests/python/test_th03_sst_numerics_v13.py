"""V13: New SST solver choices stay off by default and survive HDF5 export."""
import json
from pathlib import Path
import h5py
import pytest
from scripts.generate_th03_openfoam_blockmesh import make_case_setup, make_structured_mesh, write_case


def test_v13_defaults_are_bitwise_opt_in():
    case = make_case_setup(2, 0.25, 0.1)
    raw = case.numerics.raw_settings
    assert raw["sst_omega_newton_sink"] == "false"
    assert raw["sst_sequential_omega_first"] == "false"
    assert raw["sst_wall_omega_factor"] == "6"


def test_v13_controls_export_to_real_case(tmp_path: Path):
    dest = tmp_path / "sst-13.h5"
    write_case(dest, make_structured_mesh("smoke"), mesh_source="V13 A/B fixture",
        max_iterations=2, sst_newton_sink=True, sst_omega_first=True,
        sst_wall_factor=60, sst_inner_iterations=12)
    with h5py.File(dest, "r") as f:
        raw = json.loads(f.attrs["case_setup_json"])["numerics"]["raw_settings"]
    assert raw["sst_omega_newton_sink"] == "true"
    assert raw["sst_sequential_omega_first"] == "true"
    assert raw["sst_wall_omega_factor"] == "60"
    assert raw["sst_transport_max_iterations"] == "12"


def test_v13_wall_policy_rejects_unsupported_values():
    with pytest.raises(ValueError, match="wall factor"):
        make_case_setup(2,0.25,0.1,sst_wall_factor=10)


def test_v15_new_controls_default_off_and_export_to_case(tmp_path: Path):
    default = make_case_setup(2, 0.25, 0.1).numerics.raw_settings
    assert default['sst_compensated_residual_gate'] == 'false'
    assert default['sst_nonlinear_line_search'] == 'false'
    outfile = tmp_path / 'v15.h5'
    write_case(outfile, make_structured_mesh('smoke'), mesh_source='V15 test',
        max_iterations=2, sst_compensated_gate=True,
        sst_nonlinear_line_search=True)
    with h5py.File(outfile, 'r') as f:
        raw = json.loads(f.attrs['case_setup_json'])['numerics']['raw_settings']
    assert raw['sst_compensated_residual_gate'] == 'true'
    assert raw['sst_nonlinear_line_search'] == 'true'


def test_v15_line_search_and_fixed_damping_rejected():
    with pytest.raises(ValueError, match='exclusive'):
        make_case_setup(2, 0.25, 0.1, sst_nonlinear_line_search=True,
                        sst_picard_damping=0.5)


def test_v15_sst_linear_relaxation_validated_and_saved():
    assert make_case_setup(2,0.25,0.1).numerics.raw_settings['sst_transport_linear_relaxation'] == '0.7'
    assert make_case_setup(2,0.25,0.1,sst_linear_relaxation=1.0).numerics.raw_settings['sst_transport_linear_relaxation'] == '1.0'
    for bad in (float('nan'), float('inf'), 0.0, -0.1, 1.1):
        with pytest.raises(ValueError, match='linear relaxation'):
            make_case_setup(2,0.25,0.1,sst_linear_relaxation=bad)


def test_v16_molecular_wall_diffusion_is_explicit_opt_in():
    legacy = make_case_setup(2, 0.25, 0.1)
    proposed = make_case_setup(2, 0.25, 0.1, sst_molecular_wall_diffusion=True)
    assert legacy.numerics.raw_settings['sst_molecular_wall_diffusion'] == 'false'
    assert proposed.numerics.raw_settings['sst_molecular_wall_diffusion'] == 'true'