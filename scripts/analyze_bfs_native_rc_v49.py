"""Fail-closed, read-only V49 BFS native-flux forensic certificate.

The oracle is a diagnostic only. A successful evidence check NEVER grants
physical CFD qualification; the production nonlinear gates and independent
reference comparisons remain authoritative.
"""
from __future__ import annotations
import argparse
import json
import math
import pathlib

PREFIX = 'NATIVE_RC_ORACLE_V49 '
COMMIT = 'GLOBAL_NONLINEAR_COMMITTED_V13 '
FLOATS = ('fixed_point_defect', 'U', 'continuity', 'k', 'omega',
          'carrier_flux_gap_l1', 'carrier_flux_gap_linf',
          'native_flux_scale_linf')


def parse_record(line: str) -> dict[str, str]:
    record = {}
    for token in line.split()[1:]:
        if '=' not in token:
            continue
        key, value = token.split('=', 1)
        if key in record:
            raise ValueError(f'duplicate V49 key {key}')
        record[key] = value
    return record


def validate(log: str, reference: str | None = None,
             min_closed: int = 1, min_commits: int = 1) -> dict:
    lines = log.splitlines()
    records = [parse_record(line) for line in lines if line.startswith(PREFIX)]
    if not records:
        raise ValueError('missing V49 production oracle evidence')
    if not (min_closed >= 0 and min_commits >= 0):
        raise ValueError('negative evidence budgets')
    seen = set()
    closed_count = 0
    for record in records:
        for required in ('iteration', 'retry', 'closed', 'qualified',
                         'committed', 'physical_converged',
                         'full_six_field_oracle_qualified'):
            if required not in record:
                raise ValueError(f'missing V49 key {required}')
        key = (int(record['iteration']), int(record['retry']))
        if key in seen or key[0] < 1 or key[1] < 0:
            raise ValueError('duplicated or invalid nonlinear trial key')
        seen.add(key)
        for flag in ('qualified', 'committed', 'physical_converged',
                     'full_six_field_oracle_qualified'):
            if record[flag] != '0':
                raise ValueError(f'V49 diagnostic claimed forbidden {flag}=1')
        if record['closed'] not in ('0', '1'):
            raise ValueError('invalid V49 closed flag')
        if record['closed'] == '0':
            if 'reason' not in record and 'fixed_point_defect' not in record:
                raise ValueError('failure without diagnostic')
            continue
        closed_count += 1
        for field in FLOATS:
            if field not in record:
                raise ValueError(f'closed V49 oracle missing {field}')
            value = float(record[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'invalid V49 value for {field}')
        if float(record['fixed_point_defect']) > 1e-10:
            raise ValueError('V49 false fixed-point convergence certificate')
        for required in ('fixed_point_iterations', 'worst_face'):
            if required not in record or int(record[required]) < 0:
                raise ValueError(f'invalid V49 {required}')
        if int(record['fixed_point_iterations']) < 1:
            raise ValueError('V49 closed with no iteration')
        if record.get('repeated_same_state') not in ('0', '1'):
            raise ValueError('missing or invalid same-state repeat flag')
    if closed_count < min_closed:
        raise ValueError('not enough closed native flux evaluations')
    commits = [line for line in lines if line.startswith(COMMIT)]
    if len(commits) < min_commits:
        raise ValueError('not enough physical trial commits')
    if reference is not None:
        reference_commits = [line for line in reference.splitlines()
                             if line.startswith(COMMIT)]
        if len(reference_commits) < len(commits) or \
           reference_commits[:len(commits)] != commits:
            raise ValueError('production commitments changed under V49 probe')
    return {
        'status': 'DIAGNOSTIC_ONLY_PHYSICAL_UNQUALIFIED',
        'qualified': False,
        'physical_converged': False,
        'native_oracle_evaluations': len(records),
        'native_flux_fixed_point_closed': closed_count,
        'global_commits': len(commits),
        'production_commits_match_reference': reference is not None,
        'max_closed_native_k': max((float(rec['k']) for rec in records
                                    if rec['closed'] == '1'), default=None),
        'max_closed_native_omega': max((float(rec['omega']) for rec in records
                                        if rec['closed'] == '1'), default=None),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=pathlib.Path, required=True)
    parser.add_argument('--reference', type=pathlib.Path)
    parser.add_argument('--output', type=pathlib.Path, required=True)
    parser.add_argument('--min-closed', type=int, default=1)
    parser.add_argument('--min-commits', type=int, default=1)
    args = parser.parse_args()
    result = validate(args.log.read_text(),
                      args.reference.read_text() if args.reference else None,
                      args.min_closed, args.min_commits)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
