#!/usr/bin/env python3
"""CFDX V17: independent Rhie-Chow/patch and SST-wall evidence parser.

Reports diagnostics ONLY. It never emits a physical convergence certificate.
A converged local pressure solve does not establish global conservation.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


def fields(line: str) -> dict[str, str]:
    return dict(re.findall(r'([A-Za-z][\w]*)=([^\s]+)', line))


def finite(d: dict[str, str], key: str) -> float:
    if key not in d:
        raise ValueError(f'missing {key}')
    value = float(d[key])
    if not math.isfinite(value):
        raise ValueError(f'nonfinite {key}')
    return value


def audit_log(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    results: dict[int, dict] = {}
    wall: list[dict] = []
    pseudo: list[dict] = []
    cell_faces: dict[int, list[dict]] = {}
    pressure: list[dict] = []
    for line in path.read_text(errors='replace').splitlines():
        if line.startswith('V17_FLUX_AUDIT '):
            d = fields(line)
            idx = int(finite(d, 'iteration'))
            if idx in results:
                raise ValueError('duplicate flux audit iteration')
            data = {k: finite(d, k) for k in (
                'cell_l1', 'cell_linf', 'cell_worst', 'worst_div',
                'boundary_signed', 'boundary_abs', 'cell_signed',
                'topology_error', 'rhie_chow_identity_linf')}
            if d.get('qualified') != '0':
                raise ValueError('flux audit must not claim qualification')
            data['patches'] = {}
            results[idx] = data
        elif line.startswith('V17_FLUX_PATCH '):
            d = fields(line)
            idx = int(finite(d, 'iteration'))
            if idx not in results:
                raise ValueError('flux patch without parent audit')
            name = d.get('name')
            if not name or name in results[idx]['patches']:
                raise ValueError('missing or duplicate flux patch name')
            signed = finite(d, 'signed')
            absolute = finite(d, 'absolute')
            nfaces = finite(d, 'faces')
            if nfaces < 1 or absolute < abs(signed) - 1e-10*max(1, absolute):
                raise ValueError('invalid patch flux data')
            results[idx]['patches'][name] = {
                'faces': int(nfaces), 'signed': signed, 'absolute': absolute}
        elif line.startswith('V17_FLUX_FACE '):
            d=fields(line)
            idx=int(finite(d, 'iteration'))
            cell=int(finite(d, 'cell'))
            if idx not in results or cell!=int(results[idx]['cell_worst']):
                raise ValueError('flux worst-cell face not aligned with audit')
            sample={k: finite(d,k) for k in ('face','outward','owner_flux')}
            sample['patch']=d.get('patch')
            if d.get('qualified')!='0':
                raise ValueError('face flux diagnostic must not claim qualification')
            cell_faces.setdefault(idx,[]).append(sample)
        elif line.startswith('V17_PRESSURE_CORRECTION '):
            d=fields(line)
            sample={k:finite(d,k) for k in (
                'iteration','corrector','alpha_p','before_l1','after_l1',
                'before_linf','after_linf','reduction_ratio','pressure_linear_relative')}
            if d.get('qualified')!='0':
                raise ValueError('pressure-correction diagnostic cannot qualify')
            if (sample['alpha_p']<=0 or sample['before_l1']<0 or
                sample['after_l1']<0 or sample['before_linf']<0 or
                sample['after_linf']<0 or sample['pressure_linear_relative']<0):
                raise ValueError('invalid pressure-correction audit')
            ratio=sample['after_l1']/max(sample['before_l1'],1e-300)
            if abs(ratio-sample['reduction_ratio'])>1e-10*max(1,ratio):
                raise ValueError('pressure-correction ratio mismatch')
            pressure.append(sample)
        elif line.startswith('V17_PSEUDO_MOMENTUM '):
            d = fields(line)
            sample = {k: finite(d, k) for k in (
                'iteration', 'cfl', 'mass_min', 'mass_max',
                'physical_momentum_rhs_relative', 'continuity')}
            if d.get('qualified') != '0':
                raise ValueError('pseudo momentum cannot be marked qualified')
            if not 0 < sample['mass_min'] <= sample['mass_max'] or sample['cfl'] <= 0:
                raise ValueError('invalid pseudo momentum coefficient')
            pseudo.append(sample)
        elif line.startswith('SST_WALL_DISCRETE_AUDIT '):
            d = fields(line)
            w = {k: finite(d, k) for k in (
                'iteration', 'faces', 'd_min', 'd_max',
                'lateral_normal_max', 'omega_face_min', 'omega_face_max',
                'nasa_firstcell_min', 'nasa_firstcell_max',
                'owner_to_nasa60_min', 'owner_to_nasa60_max',
                'yplus_viscous_proxy_min', 'yplus_viscous_proxy_max')}
            if d.get('qualified') != '0' or w['faces'] <= 0 or w['d_min'] <= 0:
                raise ValueError('invalid/unqualified wall audit')
            wall.append(w)
    if not results:
        raise ValueError('no V17_FLUX_AUDIT: enable v17_flux_patch_audit')
    for i, d in results.items():
        if not d['patches']:
            raise ValueError(f'no patch audit in iteration {i}')
        signed = sum(p['signed'] for p in d['patches'].values())
        absolute = sum(p['absolute'] for p in d['patches'].values())
        scale = max(1, absolute)
        if abs(signed - d['boundary_signed']) > 1e-10*scale:
            raise ValueError(f'patch signed flux does not close in iteration {i}')
        if abs(absolute - d['boundary_abs']) > 1e-10*scale:
            raise ValueError(f'patch absolute flux does not close in iteration {i}')
        if abs(d['cell_signed'] - d['boundary_signed'] - d['topology_error']) > 1e-10*scale:
            raise ValueError(f'cell/patch flux topology mismatch in iteration {i}')
        if abs(d['topology_error']) > 1e-10*scale:
            raise ValueError(f'flux topology failure in iteration {i}')
    for i, samples in cell_faces.items():
        summed = sum(f['outward'] for f in samples)
        if abs(summed-results[i]['worst_div'])>1e-10*max(1,abs(summed)):
            raise ValueError('worst cell face fluxes do not close')
        results[i]['worst_cell_faces'] = samples
    last = results[max(results)]
    max_yplus_proxy = max((x['yplus_viscous_proxy_max'] for x in wall), default=None)
    return {
        'status': 'DIAGNOSTIC_ONLY_NOT_QUALIFIED',
        'final_iteration': max(results),
        'n_flux_audits': len(results),
        'final_flux': last,
        'max_abs_topology_error': max(abs(r['topology_error']) for r in results.values()),
        'max_rhie_chow_reconstruction_error': max(r['rhie_chow_identity_linf'] for r in results.values()),
        'pseudo': pseudo,
        'pressure_corrections': pressure,
        'wall': {'audits': len(wall), 'last': wall[-1] if wall else None,
                 'max_yplus_velocity_proxy': max_yplus_proxy,
                 'physical_yplus_qualified': False,
                 'note': 'velocity-based proxy is NOT wall-shear-derived yplus'},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--log', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    output = audit_log(args.log)
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ('pseudo',)}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())