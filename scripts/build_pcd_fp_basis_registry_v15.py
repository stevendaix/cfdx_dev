#!/usr/bin/env python3
"""Build a fail-closed coefficient-keyed Fp nullspace registry from CFP01 captures.

The basis is a PRECOMPUTED research input to the C++ PCD preconditioner.
No solver tolerance, full coupled matrix, or physical acceptance rule is changed.
For large Fp captures, this bounded dense-SVD backend refuses the request;
no guessed modes or stale basis files are used. A sparse certified backend
would be required for production-scale matrices.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import struct
import tempfile
from pathlib import Path

import numpy as np

from certify_pcd_fp_bi_nullspace_v13 import CertificateError, read_capture


def fnv64_matrix(matrix) -> int:
    n = matrix.shape[0]
    payload = (struct.pack('<QQ', n, matrix.nnz)
               + matrix.indptr.astype('<u4').tobytes()
               + matrix.indices.astype('<u4').tobytes()
               + matrix.data.astype('<f8').tobytes())
    value = 14695981039346656037
    for byte in payload:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return value


def certify_and_encode(matrix, *, rank_rtol=1e-12, min_gap=1e3):
    dense = matrix.toarray()
    u, singular, vt = np.linalg.svd(dense, full_matrices=False)
    if singular.size == 0 or not singular[0] > 0 or not np.isfinite(singular).all():
        raise CertificateError('zero/nonfinite Fp spectrum')
    rank = int(np.count_nonzero(singular > rank_rtol * singular[0]))
    nullity = len(singular) - rank
    if not (rank > 0 and 1 <= nullity <= 8):
        raise CertificateError('Fp nullspace dimension not in supported certified range 1..8')
    gap = float(singular[rank - 1] / max(singular[rank], np.finfo(float).tiny))
    left = u[:, rank:]
    right = vt[rank:, :].T
    fro = np.linalg.norm(dense)
    left_error = float(np.linalg.norm(dense.T @ left) / fro)
    right_error = float(np.linalg.norm(dense @ right) / fro)
    if not (gap >= min_gap and left_error <= 1e-12 and right_error <= 1e-12):
        raise CertificateError('Fp singular directions insufficiently separated/certified')
    fingerprint = fnv64_matrix(matrix)
    header = b'CFDXNB13' + struct.pack('<QQQ', dense.shape[0], nullity, fingerprint)
    content = header
    for basis in (left, right):
        for col in range(nullity):
            content += np.asarray(basis[:, col], dtype='<f8').tobytes()
    return content, dict(hash=f'{fingerprint:016x}', modes=nullity,
                         spectral_gap=gap, left_defect=left_error, right_defect=right_error)


def build_registry(captures: list[Path], directory: Path, max_rows: int = 1024) -> dict:
    if not captures:
        raise CertificateError('at least one Fp capture is required')
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    for capture in captures:
        # Every capture is independently checked, never treated as a reusable matrix
        # merely because it has the same nonzero pattern.
        matrix, _, _, source_hash = read_capture(capture, max_rows=max_rows)
        payload, certificate = certify_and_encode(matrix)
        destination = directory / (certificate['hash'] + '.cfpb')
        if destination.exists() and destination.read_bytes() != payload:
            # SVD signs may differ; do not silently overwrite an existing
            # certificate. A human can explicitly rebuild a fresh directory.
            raise CertificateError('existing numeric-matrix key has different basis bytes')
        if not destination.exists():
            with tempfile.NamedTemporaryFile(dir=directory, prefix='.cfpb-', delete=False) as fh:
                temporary = Path(fh.name)
                fh.write(payload)
            try:
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)
        records.append({'source': str(capture), 'capture_sha256': source_hash,
                        'basis_sha256': hashlib.sha256(payload).hexdigest(), **certificate})
    report = {'schema': 'cfdx-fp-basis-registry/v15',
              'selection': 'exact coefficient FNV64, C++ rechecks entire CSR and both nullspaces',
              'stale_fallback': False, 'physical_convergence': False,
              'dense_svd_max_rows': max_rows, 'records': records}
    (directory / 'manifest.json').write_text(
        json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capture', type=Path, action='append', required=True)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--max-rows', type=int, default=1024)
    args = parser.parse_args()
    try:
        report = build_registry(args.capture, args.directory, args.max_rows)
        print(f"PCD_FP_REGISTRY_V15 bases={len(report['records'])} stale_fallback=0 qualified=0")
        return 0
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f'PCD_FP_REGISTRY_V15_ERROR: {exc}')
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
