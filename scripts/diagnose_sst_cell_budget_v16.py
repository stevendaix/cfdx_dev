#!/usr/bin/env python3
"""Audit *signed* SST row balances on returned iterates (not a convergence certificate).

Rejects incomplete/non-finite samples and checks that convection + diffusion
- explicit source - implicit source reproduces the production A*x-b row.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

TAG = 'SST_CELL_BUDGET '
_FIELDS = ('convection', 'diffusion', 'production', 'destruction', 'sum',
           'residual', 'decomposition_error', 'source_identity_error',
           'physical_production', 'physical_destruction', 'physical_cross',
           'local_relative', 'div_phi', 'wall_orth_flux',
           'wall_owner_flux', 'wall_molecular_flux')


def inspect(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    samples: list[dict] = []
    for line in path.read_text(errors='replace').splitlines():
        if not line.startswith(TAG):
            continue
        pairs = dict(re.findall(r'([a-zA-Z_]+)=([^\s]+)', line[len(TAG):]))
        for required in ('iteration', 'field', 'cell', 'wall_faces', 'qualified', *_FIELDS):
            if required not in pairs:
                raise ValueError(f'missing {required} in SST cell budget')
        if pairs['field'] not in ('k', 'omega') or pairs['qualified'] != '0':
            raise ValueError('unexpected turbulence field or physical qualification')
        sample = {'iteration': int(pairs['iteration']), 'field': pairs['field'],
                  'cell': int(pairs['cell']), 'wall_faces': int(pairs['wall_faces'])}
        sample.update({k: float(pairs[k]) for k in _FIELDS})
        if any(not math.isfinite(sample[k]) for k in _FIELDS):
            raise ValueError('non-finite SST discrete row budget')
        source_sum = sum(sample[k] for k in ('convection', 'diffusion', 'production', 'destruction'))
        scale = max(1.0, *(abs(sample[k]) for k in _FIELDS))
        if abs(source_sum-sample['sum']) > 2e-10*scale:
            raise ValueError('invalid SST budget decomposition sum')
        if abs(sample['residual']-sample['sum']) > 2e-10*scale:
            raise ValueError('production matrix and row budget disagree')
        source_scale=max(1.0,abs(sample['physical_production']),
                         abs(sample['physical_destruction']),abs(sample['physical_cross']))
        if abs(sample['source_identity_error']) > 2e-10*source_scale:
            raise ValueError('SST frozen source differs from physical nonlinear source')
        sample['normalized_decomposition_error'] = abs(sample['decomposition_error']) / scale
        if sample['wall_faces'] > 0 and abs(sample['wall_molecular_flux']) > 1e-25:
            sample['wall_owner_molecular_ratio'] = (
                sample['wall_owner_flux'] / sample['wall_molecular_flux'])
        samples.append(sample)
    if not samples:
        raise ValueError('no SST_CELL_BUDGET samples present')
    certified = sum(line.startswith('SST_CERTIFIED_INNER')
                    for line in path.read_text(errors='replace').splitlines())
    terminal = re.findall(r'^Converged (YES|NO) Iterations (\d+)',
                          path.read_text(errors='replace'), re.M)
    return {
        'status': 'DIAGNOSTIC_ONLY_NOT_QUALIFIED',
        'sample_count': len(samples),
        'max_normalized_decomposition_error': max(s['normalized_decomposition_error']
                                                  for s in samples),
        'internal_sst_certifications': certified,
        'terminal': {'converged': terminal[-1][0] == 'YES',
                     'outer_iterations': int(terminal[-1][1])} if terminal else None,
        'samples': samples,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = inspect(args.log)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k,v in result.items() if k != 'samples'}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())