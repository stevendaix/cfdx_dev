"""Independent NASA/TMR BFS observable comparisons; never infer validation from KSP success.

Reads *raw* Driver/Seegmiller experimental files from the TMR repository.
Expected: cf.exp.dat, cp.expnew.dat and profiles.exp.dat, copied without changes
from https://tmbwg.github.io/turbmodels/Backstep_validation/ .

The result is DIAGNOSTIC_ONLY until a converged U,p,k,omega certificate,
three-grid GCI, faithful wall stresses and measured inlet Uref are proven.
Published profile data do not provide pointwise uncertainties: report errors,
never invent experimental uncertainty or an automatic physical PASS.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

SOURCE = 'https://tmbwg.github.io/turbmodels/Backstep_validation/'
STATIONS = (-4.0, 1.0, 4.0, 6.0, 10.0)

class ComparisonError(ValueError):
    pass


def _numbers(line: str, count: int) -> list[float] | None:
    try:
        parts = [float(item) for item in line.split()]
    except ValueError:
        return None
    if len(parts) != count or any(not math.isfinite(x) for x in parts):
        return None
    return parts


def read_nasa_cf(path: Path) -> list[dict]:
    rows=[]
    for line in path.read_text(encoding='utf-8').splitlines():
        values=_numbers(line,3)
        if values is None: continue
        x,cf,uncertainty=values
        if uncertainty <= 0: raise ComparisonError('invalid NASA Cf uncertainty')
        rows.append(dict(x_h=x,cf=cf,sigma=uncertainty))
    if len(rows)<10 or any(b['x_h']<=a['x_h'] for a,b in zip(rows,rows[1:])):
        raise ComparisonError('NASA Cf file truncated or nonmonotone')
    return rows


def read_nasa_cp(path: Path) -> dict:
    groups={'bottom wall':[],'top wall':[]}
    zone=None
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.strip().lower().startswith('zone'):
            match=re.search(r'"(.*?)"',line)
            zone=match.group(1).lower() if match else None
            continue
        values=_numbers(line,3)
        if values is None: continue
        if zone not in groups: raise ComparisonError('unknown NASA Cp wall patch')
        x,shifted,original=values
        groups[zone].append(dict(x_h=x,cp_shifted=shifted,cp_original=original))
    if any(len(rows)<5 for rows in groups.values()):
        raise ComparisonError('NASA Cp file incomplete')
    return groups


def read_nasa_profiles(path: Path) -> dict[float,list[dict]]:
    groups={}; zone=None
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.strip().lower().startswith('zone'):
            match=re.search(r'x/H\s*=\s*(-?[0-9.]+)',line,re.I)
            zone=float(match.group(1)) if match else None
            if zone is not None:
                if zone in groups: raise ComparisonError('duplicate NASA profile station')
                groups[zone]=[]
            continue
        vals=_numbers(line,11)
        if vals is None: continue
        if zone is None: raise ComparisonError('NASA profile row without zone')
        _,y,u,v,uu,vv,uv,*_=vals
        groups[zone].append(dict(y_h=y,u_ref=u,v_ref=v,
                                uu_ref2=uu/1000.,vv_ref2=vv/1000.,uv_ref2=uv/1000.))
    if any(x not in groups or len(groups[x])<4 for x in STATIONS):
        raise ComparisonError('NASA profile stations absent or truncated')
    return groups


def _interpolate(rows: list[dict], coordinate: str, value: str, at: float) -> float:
    if len(rows)<2 or any(not math.isfinite(float(row[k])) for row in rows for k in (coordinate,value)):
        raise ComparisonError('nonfinite or incomplete prediction samples')
    ordered=sorted(rows,key=lambda r:float(r[coordinate]))
    if any(float(b[coordinate])<=float(a[coordinate]) for a,b in zip(ordered,ordered[1:])):
        raise ComparisonError('duplicate or nonmonotone computed coordinates')
    if not float(ordered[0][coordinate])<=at<=float(ordered[-1][coordinate]):
        raise ComparisonError('computed domain does not cover experimental station')
    for a,b in zip(ordered,ordered[1:]):
        x0,x1=float(a[coordinate]),float(b[coordinate])
        if x0<=at<=x1:
            alpha=(at-x0)/(x1-x0)
            return (1-alpha)*float(a[value])+alpha*float(b[value])
    raise ComparisonError('interpolation unavailable')


def compare_cf(nasa: list[dict], predicted: list[dict]) -> dict:
    z=[]; differences=[]
    for datum in nasa:
        model=_interpolate(predicted,'x_h','cf',datum['x_h'])
        delta=model-datum['cf']; differences.append(delta)
        z.append(abs(delta)/datum['sigma'])
    return {'n':len(nasa),'max_abs_experimental_sigma':max(z),
            'rms_experimental_sigma':math.sqrt(sum(v*v for v in z)/len(z)),
            'max_absolute_error':max(abs(v) for v in differences),
            'all_points_within_3_sigma':all(v<=3 for v in z),
            'qualified':False}


def compare_cp(nasa: dict, predicted: dict, convention: str) -> dict:
    # NASA publishes two Cp origins. A good-looking curve with the wrong
    # reference pressure is NOT an experimental agreement certificate.
    if convention not in ('original','shifted'):
        raise ComparisonError('Cp convention must identify original or shifted NASA pressure origin')
    if not isinstance(predicted,dict):
        raise ComparisonError('predicted Cp must contain bottom wall and top wall patches')
    field='cp_'+convention
    comparison={}
    for wall,reference in nasa.items():
        if wall not in predicted or not isinstance(predicted[wall],list):
            raise ComparisonError(f'missing modeled Cp wall {wall}')
        differences=[_interpolate(predicted[wall],'x_h','cp',point['x_h'])-point[field]
                     for point in reference]
        comparison[wall]={'n':len(reference),
                          'max_abs_error':max(abs(delta) for delta in differences),
                          'rms_error':math.sqrt(sum(delta*delta for delta in differences)/len(differences))}
    return {'pressure_origin':convention,'walls':comparison,
            'reference_pressures_may_not_be_swapped':True,'qualified':False}


def compare_reattachment(predictions: dict) -> dict:
    value=predictions.get('reattachment_x_h')
    if not isinstance(value,(int,float)) or isinstance(value,bool) or not math.isfinite(value):
        raise ComparisonError('reattachment x/H is missing or nonfinite')
    if predictions.get('reattachment_method') != 'wall_shear_zero_crossing':
        raise ComparisonError('reattachment must come from certified wall shear sign transition')
    delta=float(value)-6.26
    return {'experimental_x_h':6.26,'experimental_sigma':0.10,
            'computed_x_h':value,'delta_x_h':delta,
            'within_experimental_band':abs(delta)<=0.10,'qualified':False}


def compare_profiles(nasa: dict[float,list[dict]], predicted: dict) -> dict:
    report={}
    for station in STATIONS:
        key=f'{station:g}'
        if key not in predicted or not isinstance(predicted[key],list):
            raise ComparisonError(f'missing modeled profile at x/H={key}')
        actual=predicted[key]
        fields={}
        for field in ('u_ref','v_ref','uv_ref2'):
            defects=[]
            for ref in nasa[station]:
                model=_interpolate(actual,'y_h',field,ref['y_h'])
                defects.append(model-ref[field])
            fields[field]={'max_abs_error':max(abs(x) for x in defects),
                           'rms_error':math.sqrt(sum(x*x for x in defects)/len(defects))}
        report[key]=fields
    return {'profiles':report,'pointwise_uncertainties_supplied':False,
            'automatic_profile_physical_pass':False}


def audit_data(nasa_dir: Path, predictions: dict | None) -> dict:
    paths=[nasa_dir/name for name in ('cf.exp.dat','cp.expnew.dat','profiles.exp.dat')]
    if not all(path.is_file() for path in paths):
        raise ComparisonError('NASA Cf/Cp/velocity-Reynolds raw files are all required')
    cf=read_nasa_cf(paths[0]); cp=read_nasa_cp(paths[1]); profiles=read_nasa_profiles(paths[2])
    output={'schema':'CFDX-BFS-NASA-RAW/v16',
            'reference_root':SOURCE,
            'raw_sha256':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in paths},
            'reattachment_exp_x_h':6.26,'reattachment_exp_sigma':0.10,
            'experimental_cf_points':len(cf),
            'experimental_cp_points':{k:len(v) for k,v in cp.items()},
            'experimental_profiles':{f'{k:g}':len(v) for k,v in profiles.items()},
            'physical_validation':'NOT_QUALIFIED'}
    if predictions is None:
        output['reason']='NO_CONVERGED_CFDX_RESULT_SUPPLIED'
        return output
    if not isinstance(predictions,dict) or any(k not in predictions for k in ('cf','cp','profiles')):
        raise ComparisonError('predicted Cf, Cp, velocity and Reynolds stress fields mandatory')
    if predictions.get('cf_method') in ('first_order','first_order_gradient',None):
        raise ComparisonError('first-order or undocumented wall shear cannot validate Cf')
    output['reattachment_diagnostic']=compare_reattachment(predictions)
    output['cp_diagnostic']=compare_cp(cp,predictions['cp'],predictions.get('cp_convention'))
    output['cf_diagnostic']=compare_cf(cf,predictions['cf'])
    output['profile_diagnostic']=compare_profiles(profiles,predictions['profiles'])
    # Even a coincidentally correct Cf cannot substitute 3-grid convergence,
    # provenance of reference velocity/Cp offset and independent nonlinear proof.
    output['reason']='NASA_PROFILE_UNCERTAINTY_AND_THREE_GRID_CERTIFICATE_STILL_REQUIRED'
    return output


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nasa-dir',type=Path,required=True)
    p.add_argument('--prediction-json',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    try:
        prediction=json.loads(args.prediction_json.read_text()) if args.prediction_json else None
        report=audit_data(args.nasa_dir,prediction)
        rc=3  # no full physical certificate is implemented or fabricated
    except (ComparisonError,OSError,ValueError,TypeError,KeyError,IndexError) as exc:
        report={'schema':'CFDX-BFS-NASA-RAW/v16','physical_validation':'NOT_QUALIFIED','reason':str(exc)}
        rc=3
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return rc

if __name__=='__main__':
    raise SystemExit(main())
