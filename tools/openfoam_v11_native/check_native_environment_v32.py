"""Preflight for a REAL OpenFOAM Foundation 11 execution.

This script never fabricates a native snapshot. It cannot certify omegaWallFunction
parity by itself. It deliberately refuses other OpenFOAM forks, absent runtime
binaries, and absent wmake; a proper solver instrumentation must be compiled and
executed separately before running qualify_openfoam_wall_v30.py.
"""
from __future__ import annotations
import argparse
import json
import os
import shutil
import sys
from pathlib import Path


def inspect_native_environment(env=None, which=None):
    env = os.environ if env is None else env
    which = shutil.which if which is None else which
    failures = []
    expected = {'WM_PROJECT': 'OpenFOAM', 'WM_PROJECT_VERSION': '11'}
    for key, wanted in expected.items():
        if env.get(key) != wanted:
            failures.append(f'{key} must equal {wanted!r}')
    root = env.get('WM_PROJECT_DIR', '')
    if not root or not Path(root).is_dir():
        failures.append('WM_PROJECT_DIR must reference an installed Foundation 11 tree')
    executables = {name: which(name) for name in ('wmake', 'foamRun')}
    for name, path in executables.items():
        if not path:
            failures.append(f'{name} executable missing')
    return {'status': 'NATIVE_RUNTIME_AVAILABLE_NOT_YET_QUALIFIED' if not failures
            else 'NATIVE_RUNTIME_UNAVAILABLE',
            'qualified': False, 'missing_requirements': failures,
            'tools': executables,
            'note': 'Native execution + same-state mesh/field snapshot and matrix parity are still required'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    result = inspect_native_environment()
    message = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(message, encoding='utf-8')
    print(message, end='')
    return 0 if result['status'] == 'NATIVE_RUNTIME_AVAILABLE_NOT_YET_QUALIFIED' else 3


if __name__ == '__main__':
    sys.exit(main())
