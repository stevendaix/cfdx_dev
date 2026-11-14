"""Fail-closed BFS V45 evidence: compare only COMMITTED reassembled physical residuals.

This is evidence for a bounded experimental Picard continuation, never a
certificate of a fully converged BFS/NASA case. Each log must contain the
initial A,b,x0 fingerprint and an explicit final execution marker.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
from compare_bfs_period4_v44 import extract, parse_pairs

_COMMITTED = re.compile(r"^GLOBAL_NONLINEAR_COMMITTED_V13\s+(.*)$", re.M)
_ESCAPE = re.compile(r'^COUPLED_PHASE_ESCAPE_V45\s+(.*)$', re.M)
_PERIOD = re.compile(r'^COUPLED_PERIOD4_FULL_V44\s+(.*)$', re.M)
_FINGERPRINT = re.compile(r'^COUPLED_MATRIX_FINGERPRINT_V40\s+(.*)$', re.M)


def compare(baseline: str, variant: str, *, requested_beta: float,
            steps: int, require_window: bool = True) -> dict:
    if not (0.0 < requested_beta <= 1.0 and math.isfinite(requested_beta)):
        raise ValueError('requested V45 beta must be in (0,1]')
    if not (1 <= steps <= 128):
        raise ValueError('requested V45 window outside 1..128')
    a, b = extract(baseline), extract(variant)
    if a['initial'] != b['initial']:
        raise ValueError('initial A,b,x0/shape fingerprint mismatch')
    # Verify every independent physical component, not only max(components).
    def physical_rows(text):
        rows = {}
        for match in _COMMITTED.finditer(text):
            pairs = parse_pairs(match.group(1))
            i = int(pairs['iteration'])
            vals = tuple(float(pairs[k]) for k in
                         ('U','p_continuity','k','omega','nut_closure','merit'))
            if i in rows or any(not math.isfinite(x) or x < 0 for x in vals):
                raise ValueError('duplicate/nonfinite physical component row')
            if not math.isclose(vals[-1], max(vals[:-1]), rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError('merit inconsistent with reassembled physical components')
            rows[i] = vals
        return rows
    physical_a, physical_b = physical_rows(baseline), physical_rows(variant)
    events = [parse_pairs(m.group(1)) for m in _ESCAPE.finditer(variant)]
    if require_window and not events:
        raise ValueError('no V45 carrier window was actually executed')
    triggered = [e for e in events if e.get('triggered') == '1']
    if require_window and not triggered:
        raise ValueError('V45 damping was not triggered by a full committed orbit')
    start = int(triggered[0]['iteration']) if triggered else None
    if len(triggered) > 1:
        raise ValueError('multi-trigger V45 run requires window-by-window evidence')
    if start is not None:
        if start not in b['trigger_iterations']:
            raise ValueError('V45 damping not tied to coherent V44 full-phase predictor')
        sample = b['commits']
        period = re.search(r'^OUTER_FIELD_PERIOD_V43 iteration=' +
                           str(start - 1) + r'\b([^\n]*)$', variant, re.M)
        if not period or 'period_four_candidate=1' not in period.group(1):
            raise ValueError('V45 trigger without preceding COMMITTED eight-field cycle')
        if (not all(i in a['commits'] for i in range(1,start)) or
            not all(i in sample for i in range(1,start))):
            raise ValueError('missing pre-intervention history')
        for i in range(1,start):
            old, new = physical_a[i], physical_b[i]
            if any(abs(t-r) > 1e-10 * max(1.0,abs(r))
                   for r,t in zip(old,new)):
                raise ValueError('V45 changed a physical state BEFORE full-cycle detection')
    for event in events:
        if (event.get('linear_candidate_certified') != '1' or
            event.get('updated_linear_solution_certified') != '0' or
            event.get('physical_qualified') != '0'):
            raise ValueError('V45 claimed a damped candidate was a linear/physical certificate')
        i = int(event['iteration'])
        beta = float(event['beta'])
        remaining = int(event['remaining'])
        if not math.isfinite(beta) or not 0 < beta <= requested_beta:
            raise ValueError('V45 trial beta is not a bounded nonlinear damping')
        if start is None or i < start or i >= start + steps:
            raise ValueError('V45 carrier damping leaked outside bounded window')
        if remaining != start+steps-i:
            raise ValueError('V45 remaining window contradicts accepted-state count')
        if i > start and event.get('triggered') != '0':
            raise ValueError('V45 re-trigger was not accounted for')
        if i in b['commits'] and event.get('retry') == '0' and not math.isclose(
                beta, requested_beta, rel_tol=1e-12,abs_tol=1e-15):
            raise ValueError('first-attempt V45 beta differs from requested value')
    if require_window and start not in b['commits']:
        raise ValueError('first V45 trial not physically committed')
    common = sorted(set(a['commits']) & set(b['commits']))
    if not common or common != list(range(1,max(common)+1)):
        raise ValueError('no contiguous common committed-state comparison')
    last = common[-1]
    post = [i for i in sorted(b['commits']) if start is not None and i >= start]
    if require_window and not post:
        raise ValueError('no post-trigger physical residual')
    final = b['commits'][max(b['commits'])]
    envelope = min(b['commits'][i] for i in post) if post else None
    return dict(status='INCOMPLETE_UNQUALIFIED',physical_converged=False,
                initial_A_b_x0=a['initial'],first_trigger=start,
                baseline_commits=len(a['commits']),trial_commits=len(b['commits']),
                common_last=last,baseline_merit_common_last=a['commits'][last],
                trial_merit_common_last=b['commits'][last],
                trial_final_merit=final,post_trigger_min_merit=envelope,
                window_events=len(events),requested_beta=requested_beta,
                steps=steps)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--trial',type=Path,required=True)
    p.add_argument('--beta',type=float,required=True)
    p.add_argument('--steps',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    report=compare(a.baseline.read_text(),a.trial.read_text(),
                   requested_beta=a.beta, steps=a.steps)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(f"V45_COMPARE status={report['status']} first_trigger={report['first_trigger']} "
          f"trial_commits={report['trial_commits']} merit={report['trial_final_merit']:.8g}")

if __name__ == '__main__':
    main()
