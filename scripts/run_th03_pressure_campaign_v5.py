"""Fail-closed three-grid TH-03 pressure campaign; never guesses Cf, Cp or xr/h."""
from __future__ import annotations
import argparse
import json
import re
import subprocess
from pathlib import Path

from generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
from audit_th03_case_contract_v10 import audit as audit_input_contract

_RESIDUAL = re.compile(r'pressure-correction solve did not converge .*?iterations=(\d+), residual=([\deE.+-]+), relative=([\deE.+-]+)')
_ACCEPTED = re.compile(r'INCOMPRESSIBLE_ITER iter=(\d+).*?pressure_linear_iterations=(\d+) pressure_relative_residual=([\deE.+-]+)')


def campaign(solver: Path, output: Path, timeout: float, tolerance: float,
             grids: tuple[str, ...]) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for grid in grids:
        mesh = make_structured_mesh(grid)
        case = output / f'th03_{grid}.cfdx.h5'
        write_case(case, mesh, mesh_source=f'OpenFOAM six-block reconstruction ({grid}); not exact blockMesh',
                   max_iterations=2, pressure_krylov='cg',
                   pressure_preconditioner='native_amg', pressure_tolerance=tolerance)
        input_contract = audit_input_contract(case)
        if input_contract['mesh_status'] == 'FAIL_GEOMETRY':
            raise ValueError(f'TH-03 mesh geometry invalid on {grid}')
        cmd = [str(solver), '--mesh', str(case), '--output-dir', str(output / f'out_{grid}'),
               '--iterations', '2', '--experimental-sst', '--benchmark-export']
        timed_out = False
        try:
            completed = subprocess.run(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, timeout=timeout, check=False)
            code = completed.returncode
            log = completed.stdout
        except subprocess.TimeoutExpired as ex:
            timed_out = True
            code = None
            log = ex.stdout or b''
            log = log.decode('utf-8', errors='replace') if isinstance(log, bytes) else log
        (output / f'{grid}.log').write_text(log)
        failures = _RESIDUAL.findall(log)
        accepted = _ACCEPTED.findall(log)
        results.append({
            'grid': grid, 'cells': mesh['n_cells'], 'timeout': timed_out,
            'input_contract_checks': input_contract['checks'],
            'inlet_eddy_viscosity_ratio': input_contract['physics']['inlet_nut_over_nu'],
            'experimental_input_comparability': input_contract['experimental_comparability'],
            'return_code': code,
            'last_pressure_failure': {
                'iterations': int(failures[-1][0]),
                'absolute_residual': float(failures[-1][1]),
                'relative_residual': float(failures[-1][2]),
            } if failures else None,
            'last_accepted_pressure_iteration': {
                'outer_iteration': int(accepted[-1][0]),
                'linear_iterations': int(accepted[-1][1]),
                'relative_residual': float(accepted[-1][2]),
            } if accepted else None,
            'physical_qualification': 'NOT_QUALIFIED',
            'log': f'{grid}.log',
        })
    summary = {
        'benchmark': 'TH-03',
        'input_lineage': 'OpenFOAM v2606 six-block reconstruction, not exact native blockMesh',
        'pressure_target': tolerance, 'results': results,
        'spatial_convergence': 'NOT_AVAILABLE',
        'cf_cp_reattachment': 'NOT_QUALIFIED',
        'reason': 'No accepted SST/2nd-order/multi-grid physically converged fields or matched inlet profiles',
    }
    (output / 'campaign.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--timeout', type=float, default=90.)
    parser.add_argument('--pressure-tolerance', type=float, default=1e-11)
    parser.add_argument('--grids', nargs='+', choices=['smoke','coarse','medium','full'],
                        default=['smoke','coarse','full'])
    args = parser.parse_args()
    if not (0 < args.timeout <= 86400.):
        parser.error('timeout must be positive and at most 86400 seconds')
    data = campaign(args.solver.resolve(), args.output, args.timeout,
                    args.pressure_tolerance, tuple(args.grids))
    print(json.dumps(data, indent=2))
    return 2  # No physically qualified three-grid campaign exists in this patch.

if __name__ == '__main__':
    raise SystemExit(main())
