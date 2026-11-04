"""Reproducible V7 SIMPLE/Rhie--Chow A/B tests; never certify truncated SST fields.

The pressure tolerance is strict (1e-11) unless mode=diagnostic.  A
successful executable exit is still not physical qualification: the source
mesh/inlet, resolved-wall treatment, conservative budgets and the 3-grid Cf/Cp
comparison must be proven separately.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import time
from pathlib import Path

import h5py

if __package__:
    from .generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
else:  # direct invocation from scripts/
    from generate_th03_openfoam_blockmesh import make_structured_mesh, write_case

PRESSURE_FAIL = re.compile(
    r'pressure-correction solve did not converge .*?iterations=(\d+), '
    r'residual=([\deE.+-]+), relative=([\deE.+-]+)'
)
PROVISIONAL_STEP = re.compile(r'SST_PROVISIONAL_OUTER iteration=(\d+)')
OUTER_ITERATION = re.compile(r'INCOMPRESSIBLE_ITER iter=(\d+)')
SST_FAIL = re.compile(
    r'SST transport did not satisfy convergence and finite-state gates '
    r'\(iterations=(\d+), converged=(\d+), k_residual=([\deE.+-]+), '
    r'omega_residual=([\deE.+-]+)'
)


def compare(solver: Path, out: Path, grids: tuple[str, ...], *,
            timeout: float = 45.0, diagnostic: bool = False,
            picard_damping: float = 0.2,
            outer_provisional: bool = False) -> dict:
    if not grids or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('nonempty grids and positive finite timeout required')
    if not (math.isfinite(picard_damping) and 0 < picard_damping <= 1):
        raise ValueError('invalid Picard damping')
    out.mkdir(parents=True, exist_ok=True)
    target = 1e-9 if diagnostic else 1e-11
    result = {'benchmark': 'TH-03', 'mode': 'DIAGNOSTIC_ONLY' if diagnostic else 'STRICT',
              'pressure_tolerance': target,
              'sst_provisional_outer': outer_provisional,
              'qualification': 'NOT_QUALIFIED', 'gci': 'UNAVAILABLE',
              'cf_cp_xr': 'UNAVAILABLE', 'entries': []}
    for grid in grids:
        mesh = make_structured_mesh(grid)
        for policy in ('arithmetic', 'projected'):
            tag = f'{grid}_{policy}'
            h5_path = out / f'{tag}.cfdx.h5'
            write_case(h5_path, mesh,
                       mesh_source='reconstructed six-block geometry; not exact published OpenFOAM polyMesh',
                       max_iterations=2, pressure_krylov='cg',
                       pressure_preconditioner='native_amg',
                       pressure_tolerance=target)
            with h5py.File(h5_path, 'r+') as h5:
                setup = json.loads(h5.attrs['case_setup_json'])
                raw = setup['numerics']['raw_settings']
                raw['pressure_face_interpolation'] = policy
                raw['sst_picard_damping'] = str(picard_damping if diagnostic else 1.0)
                raw['sst_outer_provisional'] = 'enabled' if outer_provisional else 'disabled'
                h5.attrs['case_setup_json'] = json.dumps(setup)
            cmd = [str(solver.resolve()), '--mesh', str(h5_path),
                   '--output-dir', str(out / f'output_{tag}'), '--iterations', '2',
                   '--experimental-sst', '--benchmark-export']
            t0 = time.perf_counter()
            timed_out = False
            try:
                finished = subprocess.run(cmd, capture_output=True, text=True,
                                          check=False, timeout=timeout)
                rc, output = finished.returncode, finished.stdout + '\n' + finished.stderr
            except subprocess.TimeoutExpired as exc:
                timed_out = True
                rc = None
                data = exc.stdout or b''
                output = data.decode(errors='replace') if isinstance(data, bytes) else data
            elapsed = time.perf_counter() - t0
            (out / f'{tag}.log').write_text(output, encoding='utf-8')
            p = PRESSURE_FAIL.findall(output)
            s = SST_FAIL.findall(output)
            result['entries'].append({
                'grid': grid, 'n_cells': mesh['n_cells'], 'policy': policy,
                'time_s': round(elapsed, 4), 'return_code': rc, 'timed_out': timed_out,
                'pressure_failure': {'iterations': int(p[-1][0]),
                                     'absolute_residual': float(p[-1][1]),
                                     'relative_residual': float(p[-1][2])} if p else None,
                'provisional_outer_steps': len(PROVISIONAL_STEP.findall(output)),
                'highest_logged_outer_iteration': max((int(x) for x in OUTER_ITERATION.findall(output)), default=0),
                'sst_failure': {'iterations': int(s[-1][0]),
                                'k_raw': float(s[-1][2]),
                                'omega_raw': float(s[-1][3])} if s else None,
                'physical_qualification': 'NOT_QUALIFIED',
                'log': f'{tag}.log',
            })
    (out / 'th03_v7.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--grids', nargs='+', default=['smoke','coarse','full'],
                   choices=['smoke','coarse','medium','full'])
    p.add_argument('--timeout', type=float, default=45.0)
    p.add_argument('--mode', choices=['strict','diagnostic'], default='strict')
    p.add_argument('--picard-damping', type=float, default=0.2)
    p.add_argument('--outer-provisional', action='store_true',
                   help='experimental finite/unclipped SST predictor; never certifies physics')
    args = p.parse_args()
    summary = compare(args.solver, args.output, tuple(args.grids),
                      timeout=args.timeout, diagnostic=args.mode=='diagnostic',
                      picard_damping=args.picard_damping,
                      outer_provisional=args.outer_provisional)
    print(json.dumps(summary, indent=2))
    return 2  # A/B pressure performance is NOT a physical TH-03 qualification.


if __name__ == '__main__':
    raise SystemExit(main())
