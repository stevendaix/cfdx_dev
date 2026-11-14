"""Fail-closed, read-only BFS V41 diagnostics from *actual* production logs.

This does not solve the equations or grant physical qualification; a logged
inner SST convergence and a globally committed iteration are distinct gates.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import pathlib
import re
from typing import Any

LOG_TAGS = ('SST_TWO_CYCLE_RESCUE_V41', 'SST_CERTIFIED_INNER',
            'SST_FP64_PRECISION_STALL_V14', 'GLOBAL_NONLINEAR_COMMITTED_V13',
            'COUPLED_MONOLITHIC_ILUT_V16', 'BFS_PHASE_V10')
PAIR = re.compile(r'([A-Za-z_][A-Za-z_0-9]*)=([^ ]+)')


def decode_number(x: str) -> Any:
    if x in ('true', 'false'):
        return x == 'true'
    try:
        if re.fullmatch(r'[+-]?\d+', x):
            return int(x)
        value = float(x)
        return value if math.isfinite(value) else x
    except (ValueError, OverflowError):
        return x


def parse_production_log(source: pathlib.Path) -> dict[str, Any]:
    if not source.is_file() or not source.stat().st_size:
        raise ValueError(f'missing or empty actual solver log: {source}')
    data = source.read_bytes()
    try:
        lines = data.decode('utf-8').splitlines()
    except UnicodeDecodeError as exc:
        raise ValueError(f'non-UTF-8 solver trace: {source}') from exc
    events: dict[str, list[dict[str, Any]]] = {tag: [] for tag in LOG_TAGS}
    for lineno, line in enumerate(lines, 1):
        for tag in LOG_TAGS:
            if line.startswith(tag + ' '):
                values = {k: decode_number(v) for k, v in PAIR.findall(line)}
                values['source_line'] = lineno
                events[tag].append(values)
                break
    first_certified = next(iter(events['SST_CERTIFIED_INNER']), None)
    first_rescue = next(iter(events['SST_TWO_CYCLE_RESCUE_V41']), None)
    committed = events['GLOBAL_NONLINEAR_COMMITTED_V13']
    stalled = events['SST_FP64_PRECISION_STALL_V14']
    for item in committed:
        # A committed nonlinear ITERATE is not a global convergence certificate.
        if 'qualified' not in item or 'merit' not in item:
            raise ValueError(f'incomplete global quality gate in {source}')
    return {
        'source': str(source), 'sha256': hashlib.sha256(data).hexdigest(),
        'bytes': len(data), 'first_sst_inner_certificate': first_certified,
        'first_two_cycle_rescue': first_rescue,
        'committed_outer_states': len(committed),
        'last_committed': committed[-1] if committed else None,
        'first_precision_stall': stalled[0] if stalled else None,
        'precision_stall_count': len(stalled),
        'global_physically_qualified': bool(committed and
             committed[-1].get('qualified') == 1 and
             committed[-1].get('converged') == 1),
        'ilut_recovery': events['COUPLED_MONOLITHIC_ILUT_V16'][:1],
        'phase_end': [e for e in events['BFS_PHASE_V10']
                      if e.get('phase') == 'SST_TRANSPORT_END'],
    }


def audit(actual: dict[str, pathlib.Path]) -> dict[str, Any]:
    if not actual or len(set(actual.values())) != len(actual):
        raise ValueError('at least one distinct physical log is required')
    runs = {name: parse_production_log(path) for name, path in actual.items()}
    return {'schema': 'CFDX_BFS_V41_FORENSICS_1',
            'qualification': 'diagnostics_only',
            'never_infer_physical_convergence_from_inner_sst': True,
            'runs': runs}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log', nargs=2, action='append', metavar=('LABEL', 'PATH'), required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    a = p.parse_args()
    if len(set(x[0] for x in a.log)) != len(a.log):
        p.error('duplicate run label')
    result = audit({name: pathlib.Path(path) for name, path in a.log})
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'V41_FORENSICS: {len(result["runs"])} logs, JSON={a.output}; does not certify BFS')


if __name__ == '__main__':
    main()
