"""Strict, read-only audit of experimental V51 candidate-state physical F(X).

Checks that V51 and uninstrumented production accept exactly the same states;
per-candidate replay is bitwise repeatable; the Rhie-Chow fixed point closes;
V50 scalar/vector parity is retained; and cross-couplings are observable.
No CFD physical validation is conferred by this parser.
"""
import argparse
import json
import math
from pathlib import Path


def records(text, prefix):
    out = []
    for line in text.splitlines():
        if line.startswith(prefix + ' '):
            rec = {}
            for field in line.split()[1:]:
                if '=' not in field:
                    continue
                k, v = field.split('=', 1)
                rec[k] = v
            out.append(rec)
    return out


def number(rec, key):
    if key not in rec:
        raise ValueError(f'missing {key}')
    value = float(rec[key])
    if not math.isfinite(value):
        raise ValueError(f'nonfinite {key}')
    return value


def audit(observed, baseline, min_commits=3, require_cross=False):
    base_commits = [line for line in baseline.splitlines()
                    if line.startswith('GLOBAL_NONLINEAR_COMMITTED_V13 ')]
    obs_commits = [line for line in observed.splitlines()
                   if line.startswith('GLOBAL_NONLINEAR_COMMITTED_V13 ')]
    if len(base_commits) < min_commits or obs_commits != base_commits:
        raise ValueError('V51 changed the accepted nonlinear trajectory')
    if 'FPHYS_CANDIDATE_V51 ' in baseline:
        raise ValueError('baseline contains V51 observations')
    samples = records(observed, 'FPHYS_CANDIDATE_V51')
    if len(samples) < len(obs_commits):
        raise ValueError('missing V51 candidate evaluations')
    for sample in samples:
        if sample.get('closed') != '1' or sample.get('replay_identical') != '1':
            raise ValueError('unclosed or nondeterministic V51 candidate')
        n = int(number(sample, 'cells'))
        if n <= 0 or int(number(sample, 'entries')) != 6 * n:
            raise ValueError('wrong six-field residual length')
        if not 0 <= number(sample, 'flux_defect') <= 1e-10:
            raise ValueError('V51 Rhie-Chow flux not closed at strict tolerance')
        for key in ('merit', 'continuity', 'k', 'omega', 'raw_u', 'raw_k', 'raw_omega'):
            if number(sample, key) < 0:
                raise ValueError('negative V51 residual metric')
        if sample.get('qualified') != '0' or sample.get('committed') != '0':
            raise ValueError('V51 diagnostic falsely promotes an unaccepted state')
    parity = records(observed, 'FPHYS_PARITY_V51')
    for rec in parity:
        if number(rec, 'max_normalized_gap') > 1e-9:
            raise ValueError('V51 and V50 same-state operators disagree')
        if rec.get('qualified') != '0':
            raise ValueError('false V51 parity qualification')
    cross = records(observed, 'FPHYS_CROSS_V51')
    if require_cross:
        matrix = {int(number(rec, 'column')): rec for rec in cross}
        if set(matrix) != set(range(6)):
            raise ValueError('V51 missing cross-field perturbation')
        if number(matrix[0], 'k') <= 0 or number(matrix[0], 'omega') <= 0:
            raise ValueError('V51 velocity perturbation did not alter SST')
        if number(matrix[4], 'Ux') <= 0 or number(matrix[5], 'Ux') <= 0:
            raise ValueError('V51 turbulence perturbation did not alter momentum')
        if number(matrix[3], 'continuity') <= 0:
            raise ValueError('V51 pressure perturbation did not alter continuity')
    for rec in cross:
        for key in ('Ux', 'Uy', 'Uz', 'continuity', 'k', 'omega'):
            if number(rec, key) < 0:
                raise ValueError('V51 invalid directional derivative')
        if rec.get('qualified') != '0':
            raise ValueError('false V51 directional qualification')
    return dict(status='DIAGNOSTIC_ONLY_NOT_PHYSICALLY_QUALIFIED',
                qualified=False, committed=len(obs_commits),
                attempted_candidates=len(samples), parity_samples=len(parity),
                cross_columns=len(cross), max_rhie_chow_defect=max(
                    number(s, 'flux_defect') for s in samples))


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument('--observed', type=Path, required=True)
    cli.add_argument('--baseline', type=Path, required=True)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--min-commits', type=int, default=3)
    cli.add_argument('--require-cross', action='store_true')
    options = cli.parse_args()
    result = audit(options.observed.read_text(), options.baseline.read_text(),
                   options.min_commits, options.require_cross)
    options.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
