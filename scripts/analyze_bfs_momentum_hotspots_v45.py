"""Fail-closed post-SST momentum row audit for committed BFS production states.

The solver emits each row from the actual reassembled FV operator before the
transaction commits. This tool matches those tentative rows to the subsequent
COMMITTED_V13 record, never to an unaccepted trial. 2-D BFS has active Ux, Uy;
Uz is intentionally excluded from merit even if its tiny RHS has a large
relative-to-tiny-scale value. No physical convergence is certified here.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re

_TAG = re.compile(r'^(POST_SST_MOMENTUM_ROW_V45|GLOBAL_NONLINEAR_COMMITTED_V13)\s+(.*)$')
_FINAL = re.compile(r'^Converged (YES|NO) Iterations ([0-9]+)$', re.M)
_FLOATS = ('original_signed', 'original_abs', 'rhs_scale', 'rhs_relative',
           'diagonal_term', 'off_diagonal_terms', 'rhs', 'compensated_signed',
           'pressure_source_integral', 'body_source_integral',
           'pressure_body_rhs_remainder', 'nu_t', 'k', 'omega')


def pairs(part: str) -> dict:
    result = {}
    for token in part.split():
        if '=' not in token:
            raise ValueError(f'invalid log token {token!r}')
        key,value = token.split('=', 1)
        if key in result:
            raise ValueError(f'duplicate key {key}')
        result[key] = value
    return result


def analyze(log: str, *, active_components=(0, 1), sample=(58, 59, 60, 61, 62, 63, 64, 65)) -> dict:
    if len(set(active_components)) != len(active_components) or any(i not in (0, 1, 2) for i in active_components):
        raise ValueError('invalid active momentum components')
    final = list(_FINAL.finditer(log))
    if len(final) != 1 or final[0].group(1) != 'NO':
        raise ValueError('missing or falsely successful final physical status')
    last = int(final[0].group(2))
    if last < 1:
        raise ValueError('invalid iteration count')
    pending = {}
    audited = {}
    maximum_roundoff = 0.0
    for line in log.splitlines():
        match = _TAG.match(line)
        if not match:
            continue
        kind, payload = match.groups()
        data = pairs(payload)
        iteration = int(data['iteration'])
        if not 1 <= iteration <= last:
            raise ValueError('row outside final iteration budget')
        if kind == 'POST_SST_MOMENTUM_ROW_V45':
            component = int(data['component'])
            cell = int(data['cell'])
            if not 0 <= component < 3 or cell < 0:
                raise ValueError('invalid row component or cell')
            if data.get('committed') != '0' or data.get('qualified') != '0':
                raise ValueError('tentative row incorrectly claims certification')
            values = {key:float(data[key]) for key in _FLOATS}
            if any(not math.isfinite(x) for x in values.values()):
                raise ValueError('nonfinite momentum forensic')
            if values['rhs_scale'] <= 0 or values['k'] < 0 or values['omega'] <= 0 or values['nu_t'] < 0:
                raise ValueError('nonphysical or zero scale forensic')
            if not math.isclose(values['original_abs'], abs(values['original_signed']), rel_tol=1e-10,abs_tol=1e-14):
                raise ValueError('momentum row signed and absolute mismatch')
            if not math.isclose(values['rhs_relative'], values['original_abs']/values['rhs_scale'], rel_tol=1e-9,abs_tol=1e-14):
                raise ValueError('momentum row relative scale mismatch')
            # Reevaluate the independent long-double sum from its recorded terms.
            reconstruction = (values['diagonal_term'] + values['off_diagonal_terms'] - values['rhs'])
            if not math.isclose(reconstruction,values['compensated_signed'],rel_tol=1e-8,abs_tol=1e-12):
                raise ValueError('independent CSR residual not balanced')
            if not math.isclose(values['rhs'], sum((values['pressure_source_integral'],values['body_source_integral'],values['pressure_body_rhs_remainder'])),rel_tol=1e-8,abs_tol=1e-12):
                raise ValueError('pressure and body source breakdown mismatch')
            maximum_roundoff=max(maximum_roundoff, abs(values['compensated_signed']-values['original_signed']))
            if iteration in audited:
                raise ValueError('late row after committed iteration')
            pending.setdefault(iteration,{})[component] = dict(component=component,cell=cell,**values)
        else:
            if iteration in audited:
                raise ValueError('duplicate committed iteration')
            if data.get('qualified') != '0' or data.get('converged') != '0':
                raise ValueError('unverified physical convergence claimed')
            rows = pending.pop(iteration, None)
            if rows is None or any(c not in rows for c in active_components):
                raise ValueError('committed row lacks reassembled momentum proof')
            merit = float(data['U'])
            global_merit = float(data['merit'])
            if not (math.isfinite(merit) and math.isfinite(global_merit) and merit>=0):
                raise ValueError('invalid committed physical merit')
            selected = max((rows[c] for c in active_components),key=lambda x:x['rhs_relative'])
            if not math.isclose(merit,selected['rhs_relative'],rel_tol=1e-9,abs_tol=1e-12):
                raise ValueError('row hot spot does not reproduce COMMITTED momentum gate')
            if merit > global_merit*(1+1e-9)+1e-13:
                raise ValueError('physical momentum larger than reported global merit')
            audited[iteration] = dict(merit=global_merit,U=merit,dominant_component=selected['component'],
                                      dominant_cell=selected['cell'],original_signed=selected['original_signed'],
                                      compensated_signed=selected['compensated_signed'],
                                      rhs_relative=selected['rhs_relative'],pressure_source_integral=selected['pressure_source_integral'],
                                      body_source_integral=selected['body_source_integral'],
                                      rhs=selected['rhs'],nu_t=selected['nu_t'],
                                      k=selected['k'],omega=selected['omega'])
    if len(audited) != last or sorted(audited) != list(range(1,last+1)):
        raise ValueError('missing contiguous committed-state row evidence')
    if pending:
        raise ValueError('unresolved tentative rows')
    top = Counter((x['dominant_component'],x['dominant_cell']) for x in audited.values())
    return dict(status='INCOMPLETE_UNQUALIFIED',physical_converged=False,
                iterations=last,committed_rows=len(audited),
                max_absolute_compensated_fp64_difference=maximum_roundoff,
                hotspot_counts=[dict(component=c,cell=i,states=count)
                                for (c,i),count in top.most_common(12)],
                samples={str(n):audited[n] for n in sample if n in audited})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--active-components', choices=('xy', 'xyz'), default='xy')
    args=parser.parse_args()
    report=analyze(args.log.read_text(),active_components=(0,1) if args.active_components=='xy' else (0,1,2))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('V45_ROW_EVIDENCE status='+report['status']+' committed='+str(report['committed_rows']))

if __name__=='__main__':
    main()
