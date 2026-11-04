#!/usr/bin/env python3
"""Reproducible, fail-closed TH-03 pressure and nonlinear SST diagnostic.

Never upgrades physical qualification: a pressure direct solution only audits
one linear system, while a single SST iteration cannot certify a steady flow.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

_ROW = re.compile(r'([A-Za-z_][A-Za-z_0-9]*)=([^\s]+)')


def _parsed_rows(log: str, marker: str) -> list[dict[str, str]]:
    return [dict(_ROW.findall(line)) for line in log.splitlines()
            if line.startswith(marker + ' ')]


def _number(row: dict[str, str], field: str) -> float | None:
    try:
        value = float(row[field])
    except (KeyError, ValueError):
        return None
    return value if math.isfinite(value) else None


def analyse(log_path: Path, pressure_prefix: Path | None = None,
            pressure_target: float = 1e-11) -> dict:
    log = log_path.read_text(encoding='utf-8', errors='replace')
    pressure_rows = _parsed_rows(log, 'INCOMPRESSIBLE_ITER')
    picard_rows = _parsed_rows(log, 'SST_TRANSPORT_DIAGNOSTIC')
    provisional = _parsed_rows(log, 'SST_PROVISIONAL_OUTER')
    certified_inner = _parsed_rows(log, 'SST_CERTIFIED_INNER')
    retry_count = log.count('CFDX nonlinear iteration ')  # explicit rejected steps
    linear_failed = 'pressure-correction solve did not converge' in log
    if not (pressure_rows or picard_rows or provisional or linear_failed or pressure_prefix):
        raise ValueError('no CFDX TH-03 diagnostic evidence in log')
    status_line = next((l for l in reversed(log.splitlines()) if
                        l.startswith('Converged ') or
                        l.startswith('CFDX production solver error:')), None)
    result: dict = {
        'status': 'DIAGNOSTIC_ONLY_NOT_QUALIFIED',
        'log': str(log_path),
        'solver_terminal_message': status_line,
        'pressure': {
            'iterations_reaching_carrier': len(pressure_rows),
            'direct_pressure_gate_failure': linear_failed,
            'last_pressure_relative_residual': _number(pressure_rows[-1], 'pressure_relative_residual')
            if pressure_rows else None,
            'last_continuity_normalized': _number(pressure_rows[-1], 'continuity_normalized')
            if pressure_rows else None,
        },
        'sst': {
            'certified_inner_steps': len(certified_inner),
            'last_certified_inner': certified_inner[-1] if certified_inner else None,
            'picard_samples': len(picard_rows),
            'provisional_outer_steps': len(provisional),
            'all_provisional_steps_unqualified': bool(provisional) and all(
                _number(row, 'qualified') == 0 for row in provisional),
            'rejected_outer_attempts': retry_count,
            'last_nonlinear_sample': {},
            'last_linear_contract': {},
        },
    }
    if picard_rows:
        last = picard_rows[-1]
        result['sst']['last_nonlinear_sample'] = {
            k: _number(last, k) for k in (
                'iteration', 'k_returned_nonlinear_raw', 'omega_returned_nonlinear_raw',
                'k_local_relative', 'omega_local_relative',
                'k_local_worst_cell', 'omega_local_worst_cell',
                'floor_k', 'floor_omega', 'finite_events')
        }
        result['sst']['last_linear_contract'] = {
            'k_linear_status': last.get('k_linear_status'),
            'omega_linear_status': last.get('omega_linear_status'),
            **{key: _number(last, key) for key in (
                'k_linear_relative_before_relaxation',
                'k_linear_relative_after_relaxation',
                'omega_linear_relative_before_relaxation',
                'omega_linear_relative_after_relaxation')},
        }
    if pressure_prefix is not None:
        # Existing independent SciPy/SuperLU oracle, never a replacement for
        # the CFDX physical flux and outer nonlinear convergence certificates.
        from scripts.analyse_th03_pressure_v5 import analyse as pressure_analyse
        result['pressure']['independent_sparse_lu_certificate'] = pressure_analyse(
            pressure_prefix, pressure_target)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', required=True, type=Path)
    parser.add_argument('--pressure-dump-prefix', type=Path)
    parser.add_argument('--pressure-tolerance', type=float, default=1e-11)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = analyse(args.log, args.pressure_dump_prefix, args.pressure_tolerance)
    content = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(content, encoding='utf-8')
    print(content, end='')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())