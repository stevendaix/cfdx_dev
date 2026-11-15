#!/usr/bin/env python3
"""Audit only READ-ONLY counterfactual SST mixes on BFS152 V67.

NO solver convergence certification is emitted. Missing or malformed V67
probes are a hard failure. The same historical V65 gates remain authoritative.
"""
import argparse
import json
import math
from pathlib import Path

PREFIX = 'BFS_BLOCK_COUNTERFACTUAL_V67 '
TRIAL = 'GLOBAL_SAME_OPERATOR_TRIAL_V65 '
FULL = 'BFS_FULL_SIX_COUNTERFACTUAL_V67 '
RETRIES = (0, 4, 8, 10)
FRACTIONS = ((0., 0.), (.25, .25), (.5, .5), (.75, .75),
             (1., 0.), (1., .25), (1., .5), (.5, 1.), (0., 1.))

def fields(line):
    return dict(part.split('=', 1) for part in line.split() if '=' in part)

def finite_float(fields_, key):
    value = float(fields_[key])
    if not math.isfinite(value):
        raise ValueError(f'V67 {key} is not finite')
    return value

def audit(text, require_all=True):
    rows = []
    full_six = []
    originals = {}
    for line in text.splitlines():
        if TRIAL in line:
            d = fields(line.split(TRIAL, 1)[1])
            if int(d['iteration']) == 23:
                originals[int(d['retry'])] = d
        if FULL in line:
            d = fields(line.split(FULL, 1)[1])
            if d.get('counterfactual') != '1' or d.get('committed') != '0' or d.get('qualified') != '0' or d.get('carrier_sole_beta') != '0':
                raise ValueError('V67 full-six line must be read-only and global')
            if int(d['iteration']) != 23 or int(d['retry']) != 0:
                raise ValueError('V67 full-six sweep outside bounded diagnostic')
            alpha = finite_float(d, 'alpha')
            ratio = finite_float(d, 'merit_ratio')
            growth = finite_float(d, 'largest_growth')
            if alpha <= 0. or alpha > 1. or ratio < 0. or growth < 0.:
                raise ValueError('V67 full-six invalid input')
            joint = growth <= 1.20
            monotone = ratio <= .995
            if (int(d['joint_ok']) != int(joint) or
                    int(d['monotone_ok']) != int(monotone) or
                    int(d['original_guards_ok']) != int(joint and monotone) or
                    not 0 <= int(d['worst_block']) < 6):
                raise ValueError('V67 full-six bypassed historical gates')
            full_six.append({'alpha': alpha,'ratio': ratio,'growth': growth,
                             'would_pass_original_gates': joint and monotone})
        if PREFIX in line:
            d = fields(line.split(PREFIX, 1)[1])
            if (d.get('counterfactual') != '1' or d.get('committed') != '0'
                    or d.get('transport_solved') != '0' or d.get('qualified') != '0'):
                raise ValueError('V67 probe may not certify or commit a state')
            if int(d['iteration']) != 23:
                raise ValueError('V67 probe outside bounded iteration')
            r = int(d['retry'])
            if r not in RETRIES:
                raise ValueError('V67 unrequested retry')
            beta = finite_float(d, 'beta')
            k = finite_float(d, 'k_fraction')
            omega = finite_float(d, 'omega_fraction')
            ratio = finite_float(d, 'merit_ratio')
            growth = finite_float(d, 'largest_growth')
            reference = finite_float(d, 'armijo_reference')
            if beta <= 0. or beta > 1. or ratio < 0. or growth < 0.:
                raise ValueError('V67 out-of-domain inputs')
            if not any(abs(k-kk) < 1e-14 and abs(omega-ww)<1e-14
                       for kk, ww in FRACTIONS):
                raise ValueError('V67 unexpected mixing fraction')
            if not math.isclose(reference, 1.-.005*beta, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError('V67 invalid Armijo reference')
            joint = growth <= 1.2
            monotone = ratio <= .995
            if int(d['joint_ok']) != int(joint) or int(d['monotone_ok']) != int(monotone):
                raise ValueError('V67 original individual gates incorrectly reported')
            if int(d['original_guards_ok']) != int(joint and monotone):
                raise ValueError('V67 cannot silently bypass original acceptance')
            worst = int(d['worst_block'])
            if not 0 <= worst < 6:
                raise ValueError('V67 invalid residual block')
            rows.append({'retry': r, 'beta': beta, 'k_fraction': k,
                         'omega_fraction': omega, 'ratio': ratio,
                         'growth': growth, 'would_pass_original_gates': joint and monotone,
                         'armijo_only': joint and ratio <= reference and not monotone})
    keys = {(r['retry'], r['k_fraction'], r['omega_fraction']) for r in rows}
    if len(keys) != len(rows):
        raise ValueError('V67 duplicated case or replayed log')
    if require_all:
        expected = {(r, k, w) for r in RETRIES for k, w in FRACTIONS}
        if keys != expected:
            raise ValueError('V67 probe is incomplete')
    for row in rows:
        old = originals.get(row['retry'])
        if old is None or not math.isclose(
                row['beta'], finite_float(old, 'beta'), rel_tol=1e-12, abs_tol=1e-13):
            raise ValueError('V67 did not share V65 carrier trial')
    if require_all:
        if len(full_six) != 11 or len({x['alpha'] for x in full_six}) != 11:
            raise ValueError('V67 full-six diagnostic is incomplete')
        for power in range(11):
            if not any(math.isclose(x['alpha'],2.**(-power),rel_tol=1e-12)
                       for x in full_six):
                raise ValueError('V67 full-six missing expected step')
    if full_six and 0 in originals:
        one = next(x for x in full_six if x['alpha'] == 1.0)
        old_ratio = finite_float(originals[0], 'ratio')
        if not math.isclose(one['ratio'],old_ratio,rel_tol=1e-7,abs_tol=1e-7):
            raise ValueError('V67 full-six alpha=1 must reproduce full V65 candidate')
    historical = []
    for retry, d in sorted(originals.items()):
        beta = finite_float(d, 'beta')
        ratio = finite_float(d, 'ratio')
        growth = finite_float(d, 'largest_growth') if 'largest_growth' in d else None
        if growth is not None:
            passes_blocks = growth <= 1.2
            # This is only an illustrative beta-scaled reference, not a
            # certified Armijo test, since the native V47 SST step is not
            # scaled by the same beta as the carrier state.
            historical.append({'retry': retry,'carrier_beta': beta,
                'merit_ratio': ratio,'largest_growth': growth,
                'unchanged_gate_pass': passes_blocks and ratio<=.995,
                'beta_scaled_reference_pass': passes_blocks and ratio <= 1.-.005*beta})
    return {
        'schema': 'cfdx-bfs-v67-read-only-counterfactual-2',
        'physical_qualification': False,
        'transport_solved': False,
        'compared_against_v65_unchanged': bool(originals),
        'counterfactual_count': len(rows),
        'historical_v65_trials': historical,
        'full_six_count': len(full_six),
        'full_six_original_gates_pass': sum(r['would_pass_original_gates'] for r in full_six),
        'full_six_passes': [r for r in full_six if r['would_pass_original_gates']],
        'full_six_sweep': full_six,
        'would_pass_historical_gates': sum(r['would_pass_original_gates'] for r in rows),
        'armijo_only_count': sum(r['armijo_only'] for r in rows),
        'passing_candidates': [r for r in rows if r['would_pass_original_gates']],
        'best_joint_candidate_by_retry': {
            str(retry): min((r for r in rows if r['retry'] == retry and r['growth'] <= 1.2),
                            key=lambda r: r['ratio'], default=None)
            for retry in RETRIES},
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--log', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    report = audit(args.log.read_text())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print('V67_COUNTERFACTUAL_AUDIT_COMPLETE', report['counterfactual_count'],
          'historical_gate_passes=', report['would_pass_historical_gates'],
          'physical_qualification=0')
if __name__ == '__main__':
    main()