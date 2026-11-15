@@ -0,0 +1,140 @@
#!/usr/bin/env python3
"""Independent V79 SuperLU+double-double pressure CSR replay.

No CFD solver state is committed. Output is an oracle for diagnostic full FV
Rhie-Chow re-evaluation, not a nonlinear or physical certificate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.sparse.linalg import splu


def extended_residual(a, b, hi, lo):
    ld = np.longdouble
    xv = hi.astype(ld) + lo.astype(ld)
    bv = b.astype(ld)
    data = a.data.astype(ld)
    r = np.empty(b.size, dtype=ld)
    for i in range(b.size):
        l, u = a.indptr[i], a.indptr[i + 1]
        r[i] = bv[i] - np.sum(data[l:u] * xv[a.indices[l:u]], dtype=ld)
    b2 = np.sum(bv * bv, dtype=ld)
    if not np.isfinite(b2) or b2 <= 0.:
        raise ValueError('reference RHS must be finite and nonzero')
    relative = float(np.sqrt(np.sum(r * r, dtype=ld) / b2))
    if not np.isfinite(relative):
        raise ValueError('nonfinite compensated pressure residual')
    return r, relative


def factor_and_split(matrix_path, rhs_path, tolerance=1e-11, max_passes=6):
    if not (0.0 < tolerance < 1.0) or not np.isfinite(tolerance):
        raise ValueError('original tolerance must lie in (0,1)')
    if not (0 <= max_passes <= 12):
        raise ValueError('max_passes outside 0..12')
    if np.finfo(np.longdouble).nmant <= np.finfo(np.float64).nmant:
        raise ValueError('extended long double unavailable on this architecture')
    a = mmread(str(matrix_path)).tocsr()
    a.sort_indices()
    a.sum_duplicates()
    if a.shape[0] == 0 or a.shape[0] != a.shape[1]:
        raise ValueError('matrix must be nonempty and square')
    b = np.loadtxt(rhs_path, ndmin=1).astype(np.float64)
    if b.shape != (a.shape[0],) or not np.all(np.isfinite(a.data)) or not np.all(np.isfinite(b)):
        raise ValueError('matrix/RHS size or nonfinite input')
    diag = a.diagonal()
    if np.any(diag <= 0.) or np.any(~np.isfinite(diag)):
        raise ValueError('requires an explicitly positive-diagonal pressure operator')
    lu = splu(a.tocsc(), permc_spec='COLAMD')
    high = lu.solve(b)
    if not np.all(np.isfinite(high)):
        raise ValueError('direct solve produced nonfinite pressure')
    low = np.zeros_like(high)
    fp64_relative = float(np.linalg.norm(b - a @ high) / np.linalg.norm(b))
    history = []
    for it in range(max_passes + 1):
        r, relative = extended_residual(a, b, high, low)
        history.append({'pass': it, 'relative_extended': relative})
        if relative <= tolerance or it == max_passes:
            break
        increment = lu.solve(r.astype(np.float64))
        if not np.all(np.isfinite(increment)):
            raise ValueError('nonfinite iterative refinement direction')
        exact = high.astype(np.longdouble) + low.astype(np.longdouble) + increment.astype(np.longdouble)
        next_high = exact.astype(np.float64)
        next_low = (exact - next_high.astype(np.longdouble)).astype(np.float64)
        _, next_relative = extended_residual(a, b, next_high, next_low)
        if not next_relative < relative:
            break
        high, low = next_high, next_low
    last = history[-1]['relative_extended']
    gap = a - a.T
    asymmetry = float(abs(gap).max()) if gap.nnz else 0.0
    return {
        'a': a, 'b': b, 'high': high, 'low': low,
        'report': {
            'schema': 'CFDX-PRESSURE-CSR-ORACLE-V79/1',
            'matrix_sha256': hashlib.sha256(Path(matrix_path).read_bytes()).hexdigest(),
            'rhs_sha256': hashlib.sha256(Path(rhs_path).read_bytes()).hexdigest(),
            'n': a.shape[0], 'nnz': a.nnz,
            'matrix_symmetric_absolute_defect': asymmetry,
            'max_diagonal_to_min_diagonal': float(diag.max() / diag.min()),
            'fp64_direct_relative': fp64_relative,
            'extended_split_relative': last,
            'original_tolerance': tolerance,
            'algebraically_closed': bool(last <= tolerance),
            'physical_face_flux_certified': False,
            'native_linear_certificates': False,
            'production_commit_authorized': False,
            'history': history,
        }
    }


def run(matrix_path, rhs_path, output_dir, tolerance=1e-11, max_passes=6, plot=False):
    calc = factor_and_split(matrix_path, rhs_path, tolerance, max_passes)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    report = calc['report']
    (out / 'oracle_v79.json').write_text(json.dumps(report, indent=2) + '\n')
    with (out / 'oracle_v79.split').open('w') as f:
        f.write('CFDX_PRESSURE_SPLIT_V79\n')
        f.write(f"{report['n']} {report['nnz']}\n")
        for h, l in zip(calc['high'], calc['low']):
            f.write(f'{h:.17g} {l:.17g}\n')
    if plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8.2, 4.5))
        plt.axhline(tolerance, color='black', ls='--', lw=1, label=f'Original gate: {tolerance:g}')
        plt.scatter([-1], [report['fp64_direct_relative']], label='SuperLU / double (true CSR)', marker='s', s=50)
        plt.plot([p['pass'] for p in report['history']],
                 [p['relative_extended'] for p in report['history']],
                 marker='o', label='SuperLU + high/low (extended CSR)')
        plt.yscale('log'); plt.xlabel('Refinement passes (double is x=-1)')
        plt.ylabel('Relative L2 residual on original CSR')
        plt.title('CFDX V79 - independent pressure precision (not FV convergence)')
        plt.grid(alpha=.25); plt.legend(); plt.tight_layout()
        plt.savefig(out/'pressure_csr_precision_v79.png', dpi=180)
        plt.close()
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',required=True,type=Path)
    p.add_argument('--rhs',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--tolerance',type=float,default=1e-11)
    p.add_argument('--passes',type=int,default=6)
    p.add_argument('--plot',action='store_true')
    args=p.parse_args()
    print(json.dumps(run(args.matrix,args.rhs,args.output_dir,args.tolerance,args.passes,args.plot),indent=2))

if __name__=='__main__':
    main()
