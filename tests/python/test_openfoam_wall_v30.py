"""Fail-closed offline OpenFOAM wall-function oracle regression."""
import copy
import hashlib
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile

REPO=pathlib.Path(__file__).resolve().parents[2]
SCRIPT=REPO/'scripts'/'qualify_openfoam_wall_v30.py'
spec=importlib.util.spec_from_file_location('oracle_v30',SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def pair():
    h=hashlib.sha256(b'same geometry').hexdigest()
    k=hashlib.sha256(b'same U p k omega nu_t state').hexdigest()
    ref=dict(schema='CFDX_OPENFOAM_OMEGA_V30',mesh_sha256=h,state_sha256=k,
        provenance=dict(engine='OpenFOAM Foundation',version='11',native_export=False),
        scope={'cells':[357],'faces':[10],'matrix_rows':[357]},
        cells=[dict(cell=357,omega=7.,G=2.,nut=0.3,corner_faces=2)],
        faces=[dict(face=10,owner=357,coverage=1.,weight=1.,omega=7.,G=2.,nut=0.3)],
        matrix_rows=[dict(row=357,diagonal=4.,rhs=28.,entries=[[357,4.],[363,0.]])])
    got=copy.deepcopy(ref)
    got['provenance']={'engine':'CFDX','version':'V30'}
    return ref,got

def test_identical_synthetic_is_not_certified():
    a,b=pair()
    result=mod.compare(a,b)
    assert result['status']=='INCOMPLETE_NO_NATIVE_OPENFOAM_EXPORT'
    assert not result['qualified']

def test_detects_omega_and_matrix_mismatch():
    a,b=pair()
    b['cells'][0]['omega']=8.
    b['matrix_rows'][0]['entries'][0][1]=5.
    result=mod.compare(a,b)
    assert result['status']=='FAIL_NUMERICAL_MISMATCH'
    assert any('omega' in x['at'] for x in result['first_errors'])
    assert any('entries' in x['at'] for x in result['first_errors'])

def test_rejects_different_state_or_mesh():
    a,b=pair()
    b['state_sha256']='0'*64
    try:mod.compare(a,b)
    except ValueError as e:assert 'state_sha256' in str(e)
    else:raise AssertionError('different state incorrectly passed')

def test_cli_rejects_synthetic_and_malformed_matrix():
    a,b=pair()
    with tempfile.TemporaryDirectory() as d:
        ref=pathlib.Path(d)/'ref.json'; can=pathlib.Path(d)/'cfdx.json'
        ref.write_text(json.dumps(a));can.write_text(json.dumps(b))
        run=subprocess.run([sys.executable,str(SCRIPT),str(ref),str(can)],
                           capture_output=True,text=True)
        assert run.returncode==3
        assert 'INCOMPLETE_NO_NATIVE_OPENFOAM_EXPORT' in run.stdout
        b['matrix_rows'][0]['entries'].append([357,4.])
        can.write_text(json.dumps(b))
        run=subprocess.run([sys.executable,str(SCRIPT),str(ref),str(can)],
                           capture_output=True,text=True)
        assert run.returncode==3
        assert 'FAIL_INVALID_EVIDENCE' in run.stdout

def test_scope_rejects_truncated_equal_exports():
    a,b=pair()
    a['scope']['cells']=[357,426]
    b['scope']['cells']=[357,426]
    with tempfile.TemporaryDirectory() as d:
        ref=pathlib.Path(d)/'ref.json'
        ref.write_text(json.dumps(a))
        try:mod._load(ref)
        except ValueError as e:assert 'scope' in str(e) or 'missing' in str(e)
        else:raise AssertionError('an incomplete owner set was silently accepted')


def test_claim_of_native_openfoam_without_run_log_is_insufficient():
    a,b=pair()
    a['provenance'].update(native_export=True,executable='simpleFoam',
                           case='BFS',source_commit='sha-of-revision')
    result=mod.compare(a,b)
    assert not result['qualified']
    assert result['status']=='INCOMPLETE_NO_NATIVE_OPENFOAM_EXPORT'
    with tempfile.TemporaryDirectory() as d:
        ref=pathlib.Path(d)/'ref.json'; can=pathlib.Path(d)/'cfdx.json'
        ref.write_text(json.dumps(a));can.write_text(json.dumps(b))
        run=subprocess.run([sys.executable,str(SCRIPT),str(ref),str(can)],
                           text=True,capture_output=True)
        assert run.returncode==3

def test_weight_requires_openfoam_fractional_tolerance_rule():
    a,b=pair()
    a['faces'][0]['coverage']=0.55
    a['faces'][0]['weight']=0.5
    b['faces'][0]['coverage']=0.55
    b['faces'][0]['weight']=0.5
    with tempfile.TemporaryDirectory() as d:
        ref=pathlib.Path(d)/'ref.json';can=pathlib.Path(d)/'cand.json'
        ref.write_text(json.dumps(a));can.write_text(json.dumps(b))
        assert mod._load(ref)['faces'][0]['weight']==0.5
        b['faces'][0]['weight']=0.55
        can.write_text(json.dumps(b))
        try:mod._load(can)
        except ValueError as e:assert 'weight' in str(e)
        else:raise AssertionError('incorrect wall weight accepted')
