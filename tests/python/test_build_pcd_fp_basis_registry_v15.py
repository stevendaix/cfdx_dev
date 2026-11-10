from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from build_pcd_fp_basis_registry_v15 import build_registry, certify_and_encode, fnv64_matrix
from certify_pcd_fp_bi_nullspace_v13 import CertificateError


def matrix(scale: float) -> csr_matrix:
    return csr_matrix(np.array([[scale, -scale], [2 * scale, -2 * scale]], dtype=np.float64))


def make_capture(path: Path, a: csr_matrix) -> None:
    import struct
    a = a.tocsr()
    n = a.shape[0]
    raw = b'CFDXFP01' + struct.pack('<QQ', n, a.nnz)
    raw += a.indptr.astype('<u4').tobytes()
    raw += a.indices.astype('<u4').tobytes()
    raw += a.data.astype('<f8').tobytes()
    raw += np.array([1.0, 0.0], dtype='<f8').tobytes()
    raw += np.zeros(n, dtype='<f8').tobytes()
    path.write_bytes(raw)


def test_key_tracks_numerical_coefficients():
    a, b = matrix(1.0), matrix(2.0)
    assert np.array_equal(a.indptr, b.indptr)
    assert np.array_equal(a.indices, b.indices)
    assert fnv64_matrix(a) != fnv64_matrix(b)
    _, ca = certify_and_encode(a)
    _, cb = certify_and_encode(b)
    assert ca['hash'] != cb['hash']
    assert ca['modes'] == cb['modes'] == 1


def test_registry_refreshes_by_exact_matrix(tmp_path):
    a = tmp_path / 'a.cfp'
    b = tmp_path / 'b.cfp'
    make_capture(a, matrix(1.0))
    make_capture(b, matrix(2.0))
    report = build_registry([a, b], tmp_path / 'registry', max_rows=16)
    hashes = [r['hash'] for r in report['records']]
    assert len(set(hashes)) == 2
    for key in hashes:
        assert (tmp_path / 'registry' / f'{key}.cfpb').is_file()
    assert report['stale_fallback'] is False
    assert report['physical_convergence'] is False


def test_missing_nullspace_fails_closed(tmp_path):
    a = csr_matrix(np.eye(2))
    with pytest.raises(CertificateError):
        certify_and_encode(a)


def test_oversized_capture_rejected(tmp_path):
    a = tmp_path / 'big.cfp'
    make_capture(a, matrix(1.0))
    with pytest.raises(CertificateError):
        build_registry([a], tmp_path / 'registry', max_rows=1)
