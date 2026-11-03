"""Fail-closed TH-03 resolved-wall SST prequalification from native CFDX output.

This audits numerical wall resolution only.  A passing y+ preflight is NOT
proof of turbulence accuracy, mesh convergence, or experimental validation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


class WallResolutionError(ValueError):
    pass


def finite_positive(raw: str, name: str) -> float:
    try:
        value = float(raw)
    except (ValueError, TypeError) as exc:
        raise WallResolutionError(f"{name}: nonnumeric value") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise WallResolutionError(f"{name}: requires finite strictly positive value")
    return value


def audit(face_path: Path, cell_path: Path, *, wall_patches: tuple[str,...],
          density: float, viscosity: float, max_y_plus: float = 1.0) -> dict:
    if not wall_patches or len(wall_patches) != len(set(wall_patches)):
        raise WallResolutionError('wall patches must be distinct and nonempty')
    density = finite_positive(density, 'density')
    viscosity = finite_positive(viscosity, 'molecular kinematic viscosity')
    max_y_plus = finite_positive(max_y_plus, 'max acceptable estimated y+')
    with cell_path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not {'cell_id','nu_t_m2_s'}.issubset(reader.fieldnames):
            raise WallResolutionError('missing cell fields or turbulent viscosity')
        cells = {}
        for row in reader:
            try:
                owner = int(row['cell_id'])
                value = float(row['nu_t_m2_s'])
            except (TypeError, ValueError) as exc:
                raise WallResolutionError('invalid owner cell ID or turbulent viscosity') from exc
            if owner < 0 or owner in cells or not math.isfinite(value) or value < 0.0:
                raise WallResolutionError('invalid, negative, or duplicate turbulent viscosity')
            cells[owner] = value
    if not cells:
        raise WallResolutionError('missing cells')
    entries=[]
    covered_patches: set[str] = set()
    with face_path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        needed = {'patch','owner','wall_normal_distance_m','wall_shear_x_Pa',
                  'wall_shear_y_Pa','wall_shear_z_Pa','wall_shear_method'}
        if not reader.fieldnames or not needed.issubset(reader.fieldnames):
            raise WallResolutionError('missing native wall face diagnostics')
        for row in reader:
            if row['patch'] not in wall_patches:
                continue
            if row['wall_shear_method'] != 'first_order_one_sided_diagnostic':
                raise WallResolutionError('unqualified wall shear method')
            try:
                owner = int(row['owner'])
                tau = [float(row[k]) for k in ('wall_shear_x_Pa',
                                                 'wall_shear_y_Pa','wall_shear_z_Pa')]
            except (TypeError, ValueError) as exc:
                raise WallResolutionError('invalid owner or stress components') from exc
            if owner not in cells or any(not math.isfinite(t) for t in tau):
                raise WallResolutionError('nonfinite stress or unknown owner cell')
            d = finite_positive(row['wall_normal_distance_m'], 'wall-normal distance')
            yplus = d * math.sqrt(math.hypot(*tau[:2], tau[2]) / density) / viscosity
            if 'wall_y_plus_estimate' in row and row['wall_y_plus_estimate']:
                declared = float(row['wall_y_plus_estimate'])
                if not math.isfinite(declared) or abs(declared-yplus) > 1e-7*max(1.0,yplus):
                    raise WallResolutionError('native wall y+ differs from shear-derived y+')
            if 'owner_nu_t_over_nu' in row and row['owner_nu_t_over_nu']:
                declared = float(row['owner_nu_t_over_nu'])
                if not math.isfinite(declared) or abs(declared-cells[owner]/viscosity) > 1e-7*max(1.0,declared):
                    raise WallResolutionError('native owner nu_t/nu inconsistent with cell fields')
            covered_patches.add(row['patch'])
            entries.append({'patch':row['patch'], 'y_plus_estimate':yplus,
                            'owner_nu_t_over_nu':cells[owner]/viscosity})
    if not entries:
        raise WallResolutionError('no matching wall faces')
    if set(wall_patches) != covered_patches:
        raise WallResolutionError('missing requested wall patches: '
                                  + ', '.join(sorted(set(wall_patches)-covered_patches)))
    maximum = max(e['y_plus_estimate'] for e in entries)
    minimum = min(e['y_plus_estimate'] for e in entries)
    return {'case_id':'TH-03', 'status':'PREFLIGHT_PASS' if maximum <= max_y_plus else 'BLOCKED',
            'physical_qualification':'NOT_QUALIFIED',
            'wall_method':'first_order_one_sided_diagnostic',
            'n_faces':len(entries),'min_y_plus_estimate':minimum,
            'max_y_plus_estimate':maximum,
            'max_owner_nu_t_over_nu':max(e['owner_nu_t_over_nu'] for e in entries),
            'wall_y_plus_preflight_threshold':max_y_plus,
            'reason':('resolved-wall SST y+ preflight only; full physical qualification still required'
                      if maximum<=max_y_plus else
                      'first-cell wall resolution fails the resolved-wall SST preflight')}


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--face-csv',type=Path,required=True)
    parser.add_argument('--cell-csv',type=Path,required=True)
    parser.add_argument('--patch',action='append',required=True)
    parser.add_argument('--density',type=float,required=True)
    parser.add_argument('--kinematic-viscosity',type=float,required=True)
    parser.add_argument('--max-y-plus',type=float,default=1.0)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:
        result=audit(args.face_csv,args.cell_csv,wall_patches=tuple(args.patch),
                     density=args.density,viscosity=args.kinematic_viscosity,
                     max_y_plus=args.max_y_plus)
    except (OSError,WallResolutionError,ValueError) as exc:
        parser.exit(2,f'TH-03 wall resolution preflight ERROR: {exc}\n')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'TH-03 wall {result["status"]}: max y+ estimate='
          f'{result["max_y_plus_estimate"]:.4g}; '
          f'max owner nu_t/nu={result["max_owner_nu_t_over_nu"]:.4g}')
    return 0 if result['status']=='PREFLIGHT_PASS' else 2

if __name__=='__main__':
    raise SystemExit(main())
