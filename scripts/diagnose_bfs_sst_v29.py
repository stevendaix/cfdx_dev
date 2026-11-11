#!/usr/bin/env python3
"""Fail-closed SST V29 first-outer-attempt diagnostic.

Only scores the first fully recorded SST stagnation and matching iteration-10
proposal. Never certifies SST, outer-state convergence, or NASA validation.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path


def _kv(s: str) -> dict[str, str]:
    return dict(tok.split('=', 1) for tok in s.split() if '=' in tok)


def _finite(fields: dict, name: str) -> float:
    try:
        value = float(fields[name])
    except (ValueError, KeyError) as exc:
        raise ValueError(f'missing/bad numerical token: {name}') from exc
    if not math.isfinite(value):
        raise ValueError(f'non-finite token: {name}')
    return value


def first_attempt(text: str, *, expected_active: bool) -> dict:
    stall = next((line for line in text.splitlines()
                  if line.startswith('SST_FP64_PRECISION_STALL_V14 ')), None)
    if stall is None:
        raise ValueError('no completed SST first-attempt stagnation recorded')
    first = text.split(stall, 1)[0]
    tenth = next((line for line in first.splitlines() if
                  line.startswith('SST_JFNK_V22 iteration=10 ')), None)
    if tenth is None:
        raise ValueError('first attempt iteration 10 missing')
    t = _kv(tenth)
    s = _kv(stall)
    active = t.get('v29_active') == '1'
    if t.get('v29_active') not in {'0', '1'} or active != expected_active:
        raise ValueError('V29 opt-in contradicts selected configuration')
    before, after = _finite(t, 'old_merit'), _finite(t, 'new_merit')
    terminal = _finite(s, 'k_local_relative')
    target = _finite(s, 'target')
    if not (before > after >= 0 and 0 < target <= 1e-8 and
            terminal > target and s.get('qualified') == '0' and
            s.get('accepted') == '0'):
        raise ValueError('stagnation/physical merit contract broken')
    if t.get('qualified') != '0' or t.get('accepted') != '1':
        raise ValueError('first-attempt iteration 10 requires accepted internal step only')
    if t.get('selected') not in {
        'picard', 'jfnk', 'jfnk_block', 'jfnk_spatial',
        'jfnk_spatial_log_v25', 'jfnk_smooth_v26',
        'sst_subspace_trust_v27', 'sst_adaptive_trust_v28',
        'sst_pseudo_transient_v29'}:
        raise ValueError('unknown selected SST method')
    if t.get('v29_probe') not in {'0', '1'} or t.get('v29_available') not in {'0', '1'}:
        raise ValueError('missing V29 proposal status')
    proposed = t['v29_available'] == '1'
    if not active and (t['v29_probe'] != '0' or proposed or
                       t.get('selected') == 'sst_pseudo_transient_v29'):
        raise ValueError('inactive V29 evaluated/selected a proposal')
    if proposed:
        linf = _finite(t, 'v29_linf')
        l2 = _finite(t, 'v29_l2')
        base_l2 = _finite(t, 'v26_l2_before')
        if not (active and t['v29_probe'] == '1' and
                t.get('v29_reason') == 'dual_descending' and
                linf < before and l2 < base_l2):
            raise ValueError('V29 proposal fails physical double descent')
        shifted = _finite(t, 'v29_shifted')
        wall_skip = _finite(t, 'v29_wall_skip')
        if shifted <= 0 or wall_skip < 0:
            raise ValueError('pseudo-time missing shifts or invalid wall accounting')
    elif t.get('selected') == 'sst_pseudo_transient_v29':
        raise ValueError('unavailable V29 proposal chosen')
    cfl_before = _finite(t, 'v29_cfl_before')
    cfl_after = _finite(t, 'v29_cfl_after')
    if active and not (1 / 128 <= cfl_before <= 16 and
                       1 / 128 <= cfl_after <= 16):
        raise ValueError('V29 CFL outside bounds')
    if active and t['selected'] != 'sst_pseudo_transient_v29' and \
            abs(cfl_before - cfl_after) > 1e-12:
        raise ValueError('uncommitted V29 proposal changed CFL at iteration 10')
    return {
        'v29_active': active,
        'v29_probed': t['v29_probe'] == '1',
        'v29_eligible': proposed,
        'v29_selected': t['selected'] == 'sst_pseudo_transient_v29',
        'method_selected': t['selected'],
        'iteration_10_linf_before': before,
        'iteration_10_linf_after': after,
        'candidate_linf': _finite(t, 'v29_linf') if proposed else None,
        'candidate_l2': _finite(t, 'v29_l2') if proposed else None,
        'wall_owner_rows_skipped': int(_finite(t, 'v29_wall_skip')) if active else 0,
        'shifted_transport_rows': int(_finite(t, 'v29_shifted')) if active else 0,
        'first_stagnation_iteration': int(_finite(s, 'iteration')),
        'first_stagnation_linf': terminal,
        'physical_target': target,
        'sst_converged': False,
        'outer_state_accepted': False,
        'nasa_validated': False,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('log', type=Path)
    p.add_argument('--v29', choices=['on', 'off'], required=True)
    p.add_argument('--out', type=Path)
    args = p.parse_args()
    result = first_attempt(args.log.read_text(), expected_active=args.v29 == 'on')
    output = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.out:
        args.out.write_text(output)
    else:
        print(output)


if __name__ == '__main__':
    main()
