"""Fail-closed native execution checks. These do NOT test wall parity."""
from pathlib import Path
import importlib.util

PATH = Path(__file__).resolve().parents[2] / 'tools/openfoam_v11_native/check_native_environment_v32.py'
spec = importlib.util.spec_from_file_location('cfdx_native_preflight_v32', PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

def test_missing_environment_cannot_claim_native():
    result = module.inspect_native_environment({}, lambda _: None)
    assert result['status'] == 'NATIVE_RUNTIME_UNAVAILABLE'
    assert not result['qualified']
    assert len(result['missing_requirements']) >= 4

def test_wrong_fork_is_rejected_even_if_binaries_exist(tmp_path):
    d = {'WM_PROJECT': 'OpenFOAM', 'WM_PROJECT_VERSION': 'v2606',
         'WM_PROJECT_DIR': str(tmp_path)}
    result = module.inspect_native_environment(d, lambda name: '/fake/' + name)
    assert result['status'] == 'NATIVE_RUNTIME_UNAVAILABLE'
    assert 'WM_PROJECT_VERSION must equal' in ' '.join(result['missing_requirements'])

def test_fully_declared_environment_is_still_not_a_wall_parity_certificate(tmp_path):
    d = {'WM_PROJECT': 'OpenFOAM', 'WM_PROJECT_VERSION': '11',
         'WM_PROJECT_DIR': str(tmp_path)}
    result = module.inspect_native_environment(d, lambda name: '/fake/' + name)
    assert result['status'] == 'NATIVE_RUNTIME_AVAILABLE_NOT_YET_QUALIFIED'
    assert not result['qualified']
    assert result['missing_requirements'] == []
