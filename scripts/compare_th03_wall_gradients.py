"""Independent first-order vs quadratic-normal BFS wall shear diagnostic.

Uses two aligned cell-centre samples and the exact no-slip value at a planar
horizontal wall. No pressure fitting, wall functions, or 3D unstructured
reconstruction: it deliberately *fails* outside this limited geometry.
The two-cell fit estimates sensitivity to reconstruction order, NOT a
qualified Cf, a published reattachment measurement, or a grid study.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

try:
    from scripts.extract_th03_wall_diagnostics import calculate_wall_diagnostics
except ImportError:
    from extract_th03_wall_diagnostics import calculate_wall_diagnostics


class WallFitError(ValueError):
    pass


def checked_float(value: str, name: str) -> float:
    try:
        v=float(value)
    except (TypeError,ValueError) as exc:
        raise WallFitError(f'{name} must be a number') from exc
    if not math.isfinite(v):
        raise WallFitError(f'{name} must be finite')
    return v


def quadratic_dudn(d1: float, u1: float, d2: float, u2: float) -> float:
    """Wall derivative for U(0)=0, U(d)=a*d+b*d^2, 0<d1<d2."""
    if not (all(math.isfinite(x) for x in (d1,d2,u1,u2)) and 0<d1<d2):
        raise WallFitError('singular/nonfinite quadratic wall-normal stencil')
    result=(u1*d2*d2-u2*d1*d1)/(d1*d2*(d2-d1))
    if not math.isfinite(result):
        raise WallFitError('nonfinite quadratic wall derivative')
    return result


def crossing(points: list[tuple[float,float]]) -> float:
    ordered=sorted(points)
    roots=[]
    for (xa,ca),(xb,cb) in zip(ordered,ordered[1:]):
        if ca<0<cb:
            roots.append(xa-ca*(xb-xa)/(cb-ca))
    if len(roots)!=1:
        raise WallFitError(f'expected exactly one negative-positive reattachment, found {len(roots)}')
    return roots[0]


def compare(face_csv: Path, cell_csv: Path, *, patch: str, density: float,
            viscosity: float, bulk_speed: float, step_height: float,
            step_x: float = 0.0) -> dict:
    if not all(math.isfinite(x) for x in (density,viscosity,bulk_speed,step_height,step_x)) or \
       min(density,viscosity,abs(bulk_speed),step_height)<=0:
        raise WallFitError('invalid physical normalization or fluid viscosity')
    first_profile, first_root=calculate_wall_diagnostics(face_csv,patch=patch,
        density=density,bulk_speed=bulk_speed,step_height=step_height,step_x=step_x)
    with cell_csv.open(newline='',encoding='utf-8') as stream:
        cells=list(csv.DictReader(stream))
    required={'x_m','y_m','z_m','u_x_m_s'}
    if not cells or not required.issubset(cells[0]):
        raise WallFitError('cell CSV is missing geometry or velocity')
    samples=[]
    for cell in cells:
        samples.append(tuple(checked_float(cell[name],name) for name in
                        ('x_m','y_m','z_m','u_x_m_s')))
    with face_csv.open(newline='',encoding='utf-8') as stream:
        faces=list(csv.DictReader(stream))
    if not faces or not {'patch','x_m','y_m','z_m','area_x_m2','area_y_m2',
       'area_z_m2','area_m2','wall_normal_distance_m','wall_shear_x_Pa'}.issubset(faces[0]):
        raise WallFitError('face CSV is missing geometry or wall shear')
    data=[]
    for face in faces:
        if face['patch']!=patch:continue
        x,y,z,ax,ay,az,area,first_dist,tau1=(checked_float(face[name],name)
            for name in ('x_m','y_m','z_m','area_x_m2','area_y_m2',
                         'area_z_m2','area_m2','wall_normal_distance_m','wall_shear_x_Pa'))
        if (x-step_x)<0:continue
        if not (area>0 and ay<0 and abs(ax)/area<1e-10 and
                abs(az)/area<1e-10 and abs(ay/area+1)<1e-10):
            raise WallFitError('not a planar horizontal outward-down wall face')
        tol=1e-9*max(step_height,abs(x),1)
        candidate=sorted((cy-y,ux) for cx,cy,cz,ux in samples
                         if abs(cx-x)<=tol and abs(cz-z)<=tol and cy>y)
        if len(candidate)<2:
            raise WallFitError('wall station lacks two vertically aligned cell centres')
        (d1,u1),(d2,u2)=candidate[:2]
        if abs(d1-first_dist)>1e-8*max(1,first_dist):
            raise WallFitError('first cell does not match wall normal distance')
        mu=density*viscosity
        if abs(mu*u1/d1-tau1)>1e-7*max(1,abs(tau1)):
            raise WallFitError('face shear and cell state do not agree at first order')
        tau2=mu*quadratic_dudn(d1,u1,d2,u2)
        data.append({'x':(x-step_x)/step_height,
                     'cf_first':2*tau1/(density*bulk_speed**2),
                     'cf_quadratic':2*tau2/(density*bulk_speed**2),
                     'd1_m':d1,'d2_m':d2})
    if len(data)<3:raise WallFitError('not enough downstream wall stations')
    data.sort(key=lambda item:item['x'])
    quadratic_root=crossing([(p['x'],p['cf_quadratic']) for p in data])
    return {'status':'DIAGNOSTIC_ONLY','physical_qualification':'NOT_QUALIFIED',
            'first_order_reattachment_x_over_h':first_root,
            'quadratic_reattachment_x_over_h':quadratic_root,
            'reconstruction_spread_x_over_h':abs(first_root-quadratic_root),
            'quadratic_method':'two_cell_normal_quadratic_diagnostic',
            'n_wall_stations':len(data),'profiles':data}


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--face-csv',type=Path,required=True)
    p.add_argument('--cell-csv',type=Path,required=True)
    p.add_argument('--patch',default='downstream_wall')
    p.add_argument('--density',type=float,required=True)
    p.add_argument('--kinematic-viscosity',type=float,required=True)
    p.add_argument('--bulk-speed',type=float,required=True)
    p.add_argument('--step-height',type=float,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    try:
        result=compare(args.face_csv,args.cell_csv,patch=args.patch,density=args.density,
            viscosity=args.kinematic_viscosity,bulk_speed=args.bulk_speed,
            step_height=args.step_height)
    except (OSError,WallFitError,ValueError) as exc:
        p.exit(2,f'TH-03 quadratic wall-fit diagnostics FAIL: {exc}\n')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'TH-03 Cf first-order x_r/h={result["first_order_reattachment_x_over_h"]:.6g}; '
          f'quadratic x_r/h={result["quadratic_reattachment_x_over_h"]:.6g} '
          '(diagnostic only)')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
