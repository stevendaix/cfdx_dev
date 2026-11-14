"""V47 read-only fail-closed trial/rollback evidence (never a BFS qualification)."""
from __future__ import annotations
import argparse
import json
import math
import pathlib
import re
from typing import Any


def rows(text: str, tag: str) -> list[dict[str, str]]:
    return [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in text.splitlines()
            if line.startswith(tag + ' ')]


def audit(text: str, *, shrink: float = 0.5, expected_first_beta: float = 1.0,
          require_deep_retry: bool = False,
          require_complete_log: bool = False) -> dict[str, Any]:
    if not math.isfinite(shrink) or not 0 < shrink < 1:
        raise ValueError('invalid shrink')
    trials = rows(text, 'COUPLED_COMPOSITE_TRIAL_V47')
    if not trials:
        raise ValueError('no V47 carrier trials in production log')
    retries: dict[int, list[int]] = {}
    joint_by_iter: dict[int, list[float]] = {}
    for item in rows(text, 'GLOBAL_JOINT_V13'):
        try:
            if item.get('accepted')=='0':
                joint_by_iter.setdefault(int(item['iteration']), []).append(float(item['growth']))
        except (KeyError, ValueError) as exc:
            raise ValueError('invalid physical V13 rejection evidence') from exc
    for item in trials:
        try:
            i, r, beta = int(item['iteration']), int(item['retry']), float(item['beta'])
        except (KeyError, ValueError) as exc:
            raise ValueError('corrupt V47 trial') from exc
        if i < 1 or r < 0 or r > 12 or not math.isfinite(beta):
            raise ValueError('invalid iteration/retry/beta')
        geometric=expected_first_beta * (shrink ** r)
        kind=item.get('search_kind','geometric')
        if kind=='quadratic':
            values=joint_by_iter.get(i,[])
            if r!=3 or len(values)<3:
                raise ValueError('quadratic retry without three physical rejection samples')
            x0,x1,x2=(1.0,shrink,shrink*shrink)
            y0,y1,y2=values[:3]
            s01=(y0-y1)/(x0-x1)
            s12=(y1-y2)/(x1-x2)
            a=(s01-s12)/(x0-x2)
            b=s01-a*(x0+x1)
            if not a>0 or not math.isfinite(a):
                raise ValueError('invalid quadratic physical objective curvature')
            prediction=-b/(2*a)
            if not (x2+0.1*(x1-x2)<prediction<x1-0.1*(x1-x2)):
                raise ValueError('unsafeguarded polynomial step')
            expected=prediction
        elif kind=='geometric':
            expected=geometric
        else:
            raise ValueError('unrecognized V47 search kind')
        if abs(beta-expected)>1e-10:
            raise ValueError('V47 beta violates physical-trial step selection')
        if any(item.get(k) != v for k, v in {
            'carrier_linear_candidate_certified': '1',
            'damped_linear_solution_certified': '0',
            'flux_rederived': '1',
            'sst_transport_pending': '1',
            'physical_qualified': '0', 'committed': '0'}.items()):
            raise ValueError('V47 claimed unsupported state certification')
        retries.setdefault(i, []).append(r)
    for i, seq in retries.items():
        if seq != list(range(len(seq))):
            raise ValueError('noncontiguous V47 trial history at outer step %d' % i)
    if require_deep_retry and not (any(max(seq) >= 4 for seq in retries.values()) or
                                   any(t.get('search_kind')=='quadratic' for t in trials)):
        raise ValueError('insufficient backtracking or interpolated evidence')
    if require_complete_log and not re.search(r'^Converged NO Iterations [0-9]+\s*$',text,re.M):
        raise ValueError('missing explicit non-converged production termination')
    commits = rows(text, 'GLOBAL_NONLINEAR_COMMITTED_V13')
    v8 = rows(text, 'GLOBAL_NONLINEAR_V8')
    v14 = rows(text, 'GLOBAL_CONTRACTION_TRIAL_V14')
    joint = rows(text, 'GLOBAL_JOINT_V13')
    if not commits:
        raise ValueError('no committed physical state')
    seen: set[int] = set()
    previous_merit: float | None = None
    for item in commits:
        try:
            i = int(item['iteration'])
            metrics = [float(item[k]) for k in ('U', 'p_continuity', 'k', 'omega', 'nut_closure')]
            merit = float(item['merit'])
        except (ValueError, KeyError) as exc:
            raise ValueError('corrupt committed physical certificate') from exc
        if (i in seen or i not in retries or item.get('converged') != '0'
                or item.get('qualified') != '0'
                or not all(math.isfinite(v) and v >= 0 for v in metrics+[merit])):
            raise ValueError('duplicate/missing/nonfinite committed iteration')
        if abs(merit - max(metrics)) > max(1e-11, 1e-9 * merit):
            raise ValueError('global merit does not match physical components')
        # A commit needs positive V14 AND matching independently re-evaluated
        # V8 components, not just a plausible-looking V13 summary line.
        same_v8 = [v for v in v8 if v.get('iteration') == str(i)]
        accepted_v14 = [v for v in v14 if v.get('iteration') == str(i)]
        if len(same_v8) != 1 or len(accepted_v14) != 1:
            raise ValueError('missing/duplicate independent V8 or V14 evidence')
        if accepted_v14[0].get('allowed') != '1':
            raise ValueError('V14 rejected this purported committed state')
        for k, v in zip(('U', 'p_continuity', 'k', 'omega', 'nut_closure'), metrics):
            try: observed = float(same_v8[0][k])
            except (KeyError, ValueError) as exc:
                raise ValueError('malformed independent V8 evidence') from exc
            if not math.isfinite(observed) or abs(observed-v)>max(1e-11, 1e-9*v):
                raise ValueError('V8/V13 physical component disagreement')
        try: v14_merit = float(accepted_v14[0]['new'])
        except (KeyError, ValueError) as exc:
            raise ValueError('malformed V14 merit') from exc
        if not math.isfinite(v14_merit) or abs(v14_merit-merit)>max(1e-11,1e-9*merit):
            raise ValueError('V14/V13 merit disagreement')
        if i > 1:
            associated = [v for v in joint if v.get('iteration') == str(i)]
            if not associated or associated[-1].get('accepted') != '1':
                raise ValueError('missing V13 accepted trial for committed state')
        # V14 imposes strict 0.995 contraction, unless already converged.
        if previous_merit is not None and merit > 0.99500001 * previous_merit:
            raise ValueError('committed global merit violated V14')
        previous_merit = merit
        seen.add(i)
    if sorted(seen) != list(range(1, len(seen)+1)):
        raise ValueError('committed state sequence contains holes')
    finished = re.search(r'^Converged NO Iterations (\d+)\s*$',text,re.M)
    if require_complete_log and (finished is None or int(finished.group(1)) != len(commits)):
        raise ValueError('termination count disagrees with committed physical states')
    if 'Converged YES' in text:
        raise ValueError('unverified physical convergence claim')
    aborted = 'CFDX production solver error:' in text
    return {'status': 'ABORTED_UNQUALIFIED' if aborted else 'INCOMPLETE_UNQUALIFIED',
            'physical_converged': False,
            'committed_outer_steps': len(commits), 'total_trial_steps': len(trials),
            'quadratic_proposals': sum(t.get('search_kind')=='quadratic' for t in trials),
            'deepest_retry': max(max(s) for s in retries.values()),
            'last_committed_merit': previous_merit,
            'trial_counts': {str(i):len(s) for i,s in retries.items()},
            'note':'V14 physical merit contract; no physical BFS certificate'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--require-deep-retry', action='store_true')
    p.add_argument('--require-complete-log', action='store_true')
    args = p.parse_args()
    result = audit(args.log.read_text(), require_deep_retry=args.require_deep_retry,
                   require_complete_log=args.require_complete_log)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, sort_keys=True))

if __name__ == '__main__':
    main()
