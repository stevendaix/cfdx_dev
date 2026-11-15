#!/usr/bin/env python3
"""Fail-closed native seventh-Schur replay on SIX identical native BFS matrices.

The replay is algebraic on frozen 4N. It cannot claim a converged RANS/SST BFS.
Default target is 1e-9 for THIS diagnostic; this script does not alter any
CFDX physical or production linear tolerance. It never renames PCD as selfp.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys

METHODS = ('pcd', 'lsc', 'bfbt', 'block_local', 'simple', 'simplec')
PATTERN = re.compile(r'V63_SELFP_REPLAY n=(\d+) nnz=(\d+) schur_nnz=(\d+) '
                     r'iterations=(\d+) true_relative=([\deE+.-]+) '
                     r'claimed_relative=([\deE+.-]+) status=(\d+) physical_qualified=0')

def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def audit(replay, captures, timeout=60):
    if not replay.is_file() or not captures.is_dir():
        raise ValueError('replay binary or captures directory is missing')
    hashes = {kind: {} for kind in ('A', 'b', 'x0')}
    for method in METHODS:
        for kind in hashes:
            path = captures / f'{method}_{kind}.mtx'
            if not path.is_file():
                raise ValueError(f'missing physical operator: {path.name}')
            hashes[kind][method] = hash_file(path)
    for kind, per_method in hashes.items():
        if len(set(per_method.values())) != 1:
            raise ValueError(f'six-method input mismatch: {kind}; refuse to compare')
    results = {}
    for method in METHODS:
        try:
            proc = subprocess.run([str(replay),
                str(captures / f'{method}_A.mtx'),
                str(captures / f'{method}_b.mtx'),
                str(captures / f'{method}_x0.mtx')],
                capture_output=True, text=True, timeout=timeout, check=False)
        except subprocess.TimeoutExpired as exc:
            raise ValueError(f'{method}: replay timeout, no certificate') from exc
        match = PATTERN.search(proc.stdout)
        if proc.returncode or not match:
            raise ValueError(f'{method}: replay failure {proc.returncode} {proc.stderr}')
        n, nnz, snnz, count, true, claimed, status = match.groups()
        true, claimed = float(true), float(claimed)
        if int(n) != 608 or int(nnz) != 7080 or int(snnz) == 0 or \
           not (0 <= int(count) <= 500) or int(status) != 0 or \
           not math.isfinite(true) or not math.isfinite(claimed) or \
           true > 1e-9 or abs(true-claimed) > 1e-12:
            raise ValueError(f'{method}: not certified on the original 4N matrix')
        results[method] = {'iterations': int(count), 'true_relative': true,
                           'status': int(status), 'algebraic_only': True}
    return {'format': 'CFDX-NATIVE-SELFP-V63-ALGEBRAIC-ONLY',
            'same_original_A_b_x0_sha256': True,
            'source_sha256': {kind: next(iter(per.values())) for kind, per in hashes.items()},
            'results': results, 'nonlinear_BFS_qualified': False,
            'PCD_changed': False, 'six_original_schur_changed': False}

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replay', type=Path, required=True)
    p.add_argument('--captures', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(argv)
    try:
        evidence = audit(args.replay, args.captures)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(evidence,indent=2,sort_keys=True)+'\n')
        print('V63_SELFP_SIX_COMMON_INPUT_PASS all six original matrices unchanged; algebraic only')
        return 0
    except (ValueError, OSError) as e:
        print(f'V63_SELFP_SIX_COMMON_INPUT_REJECT {e}',file=sys.stderr)
        return 2

if __name__ == '__main__':
    sys.exit(main())