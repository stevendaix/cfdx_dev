"""V31 OpenFOAM-native hook packaging tests; NOT an OpenFOAM parity test.

These fixtures are SYNTHETIC BY DESIGN and must never count as native evidence.
"""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PKG_SCRIPT = ROOT / 'tools/openfoam_v11_native/package_native_snapshot_v31.py'
ORACLE_SCRIPT = ROOT / 'scripts/qualify_openfoam_wall_v30.py'

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

pkg=load(PKG_SCRIPT,'of_native_packager_v31')
oracle=load(ORACLE_SCRIPT,'of_native_oracle_v30')

def fixtures(tmp_path):
    paths={}
    for tag, content in {'mesh':'canonical mesh bytes',
                         'state':'same k omega G nut U p',
                         'run':'simulated-not-OpenFOAM-runtime-log'}.items():
        paths[tag]=tmp_path/(tag+'.dat')
        paths[tag].write_text(content)
    raw={'schema':'CFDX_OF11_NATIVE_SNAPSHOT_V31',
         'stage':'post_manipulate_pre_solve',
         'cells':[{'cell':357,'omega':7.0,'G':2.0,'nut':0.3,'corner_faces':2}],
         'faces':[{'face':10,'owner':357,'coverage':0.55,'weight':0.5,
                   'omega':7.,'G':2.,'nut':0.3}],
         'matrix_rows':[{'row':357,'diagonal':4.,'rhs':28.,
                         'entries':[[357,4.],[363,0.]]}]}
    paths['raw']=tmp_path/'raw.json'
    paths['raw'].write_text(json.dumps(raw))
    return paths,raw

def package(paths):
    return pkg.package(paths['raw'],paths['mesh'],paths['state'],paths['run'],
                       'UNVERIFIED-synthetic-fixture','BFS', 'unverified-test')

def test_packaged_fixture_is_only_schema_evidence(tmp_path):
    paths,_=fixtures(tmp_path)
    record=package(paths)
    output=tmp_path/'native_record.json'; output.write_text(json.dumps(record))
    validated=oracle._load(output)
    assert validated['scope']=={'cells':[357],'faces':[10],'matrix_rows':[357]}
    assert validated['provenance']['run_log_sha256']==hashlib.sha256(paths['run'].read_bytes()).hexdigest()
    candidate=copy.deepcopy(validated)
    candidate['provenance']={'engine':'CFDX','version':'V31'}
    assert not oracle.compare(validated,candidate)['qualified']

def test_rejects_missing_or_wrong_stage(tmp_path):
    paths,raw=fixtures(tmp_path)
    raw['stage']='before_assembly'
    paths['raw'].write_text(json.dumps(raw))
    with pytest.raises(ValueError,match='native post-manipulation'):
        package(paths)

def test_rejects_duplicate_owner_face_or_row(tmp_path):
    paths,raw=fixtures(tmp_path)
    for group in ('cells','faces','matrix_rows'):
        invalid=copy.deepcopy(raw)
        invalid[group].append(copy.deepcopy(invalid[group][0]))
        paths['raw'].write_text(json.dumps(invalid))
        with pytest.raises(ValueError,match='invalid'):
            package(paths)

def test_rejects_incomplete_rows_and_inconsistent_weight(tmp_path):
    paths,raw=fixtures(tmp_path)
    invalid=copy.deepcopy(raw)
    invalid['matrix_rows'][0].pop('entries')
    paths['raw'].write_text(json.dumps(invalid))
    record=package(paths)
    output=tmp_path/'incomplete.json'; output.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='matrix'):
        oracle._load(output)
    invalid=copy.deepcopy(raw)
    invalid['faces'][0]['weight']=0.55
    paths['raw'].write_text(json.dumps(invalid))
    record=package(paths); output.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='weight'):
        oracle._load(output)

def test_rejects_missing_native_log(tmp_path):
    paths,_=fixtures(tmp_path)
    paths['run'].unlink()
    with pytest.raises(FileNotFoundError):
        package(paths)

def test_rejects_nonfinite_and_out_of_scope_rows(tmp_path):
    paths,raw=fixtures(tmp_path)
    invalid=copy.deepcopy(raw)
    invalid['cells'][0]['omega']=float('nan')
    paths['raw'].write_text(json.dumps(invalid))
    with pytest.raises(ValueError,match='nonfinite'):
        package(paths)
    invalid=copy.deepcopy(raw)
    invalid['faces'][0]['owner']=200
    paths['raw'].write_text(json.dumps(invalid))
    record=package(paths)
    out=tmp_path/'wrong_owner.json'; out.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='omitted'):
        oracle._load(out)
