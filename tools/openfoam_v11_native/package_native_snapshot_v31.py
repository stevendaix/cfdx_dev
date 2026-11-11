"""Package a *native* OpenFOAM Foundation v11 instrumentation snapshot for
CFDX V30 oracle. Requires byte-identical canonical mesh/state reference
files and the actual native execution log. Not an OpenFOAM emulator.
"""
import argparse,hashlib,json,math
from pathlib import Path

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def package(raw,mesh_file,state_file,run_log,executable,case,commit):
    data=json.loads(Path(raw).read_text(encoding='utf-8'))
    if data.get('schema')!='CFDX_OF11_NATIVE_SNAPSHOT_V31' or data.get('stage')!='post_manipulate_pre_solve':
        raise ValueError('not a native post-manipulation OpenFOAM snapshot')
    groups=('cells','faces','matrix_rows')
    keys=('cell','face','row')
    scope={}
    for group,key in zip(groups,keys):
        entries=data.get(group)
        if not isinstance(entries,list) or not entries:
            raise ValueError(f'missing {group}')
        indices=[e.get(key) for e in entries]
        if any(type(i) is not int or i<0 for i in indices) or len(set(indices))!=len(indices):
            raise ValueError(f'invalid {group} identifiers')
        scope[group]=sorted(indices)
        for e in entries:
            for k,v in e.items():
                if isinstance(v,float) and not math.isfinite(v):
                    raise ValueError(f'nonfinite {group}.{k}')
    if not all([executable,case,commit]):
        raise ValueError('native executable, case and source commit are mandatory')
    out={'schema':'CFDX_OPENFOAM_OMEGA_V30','mesh_sha256':digest(mesh_file),
         'state_sha256':digest(state_file),'scope':scope,
         'provenance':{'engine':'OpenFOAM Foundation','version':'11',
          'native_export':True,'executable':executable,'case':case,
          'source_commit':commit,'run_log_sha256':digest(run_log)}}
    out.update({k:data[k] for k in groups})
    # Reuse the exact independent validator: not even a self-attested raw
    # snapshot is promoted when its native execution log is absent.
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
    from qualify_openfoam_wall_v30 import _load
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native_snapshot',type=Path)
    p.add_argument('--canonical-mesh',required=True,type=Path)
    p.add_argument('--canonical-state',required=True,type=Path)
    p.add_argument('--native-run-log',required=True,type=Path)
    p.add_argument('--executable',required=True)
    p.add_argument('--case',required=True)
    p.add_argument('--source-commit',required=True)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    record=package(a.native_snapshot,a.canonical_mesh,a.canonical_state,
                   a.native_run_log,a.executable,a.case,a.source_commit)
    a.output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    # Validate the output with the unchanged fail-closed V30 schema checks.
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
    from qualify_openfoam_wall_v30 import _load
    _load(a.output)
    print('Native snapshot schema and log provenance recorded; no CFDX numerical parity claimed')
if __name__=='__main__': main()
