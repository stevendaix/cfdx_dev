"""V13 numerical SVD truth contract; no PCD / BFS convergence inferred."""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from certify_pcd_fp_bi_nullspace_v13 import CertificateError, analyze, main


def fixture(tmp_path, matrix, rhs):
    arr = np.asarray(matrix, dtype="<f8")
    rhs = np.asarray(rhs, dtype="<f8")
    rows = [0]
    cols = []
    values = []
    for row in arr:
        for j, val in enumerate(row):
            if val:
                cols.append(j)
                values.append(val)
        rows.append(len(cols))
    n = len(rhs)
    raw = b"CFDXFP01" + struct.pack("<QQ", n, len(cols))
    raw += np.asarray(rows, dtype="<u4").tobytes()
    raw += np.asarray(cols, dtype="<u4").tobytes()
    raw += np.asarray(values, dtype="<f8").tobytes()
    raw += rhs.tobytes()
    raw += np.zeros(n, dtype="<f8").tobytes()
    path = tmp_path / "case.cfp"
    path.write_bytes(raw)
    return path


def test_two_sided_nonsymmetric_incompatible(tmp_path):
    path = fixture(tmp_path, [[1, -1], [2, -2]], [1, 0])
    c = analyze(path)
    assert c["matrix"]["numerical_nullity_left"] == 1
    assert c["matrix"]["numerical_nullity_right"] == 1
    assert c["numerical_basis_certified"]
    assert c["incompatible_rhs_relative"] == pytest.approx(2 / (5 ** 0.5))
    assert c["original_true_residual_relative_rhs"] == pytest.approx(c["incompatible_rhs_relative"])
    assert c["projected_true_residual_l2"] < 1e-12
    assert not c["original_rhs_compatible_at_1e-8"]
    assert main(["--capture", str(path), "--output", str(tmp_path / "proof.json"), "--require-compatible"]) == 3
    assert json.loads((tmp_path / "proof.json").read_text())["physical_convergence"] is False


def test_exact_null_compatible_rhs(tmp_path):
    path = fixture(tmp_path, [[1, -1], [2, -2]], [1, 2])
    c = analyze(path)
    assert c["original_rhs_compatible_at_1e-8"]
    assert c["original_true_residual_relative_rhs"] < 1e-12
    assert main(["--capture", str(path), "--output", str(tmp_path / "compatible.json"), "--require-compatible"]) == 0


def test_nonsingular_and_empty_nullspace(tmp_path):
    path = fixture(tmp_path, [[2, -1], [-1, 2]], [1, 0])
    c = analyze(path)
    assert c["matrix"]["numerical_nullity_left"] == 0
    assert c["original_true_residual_relative_rhs"] < 1e-12
    assert c["original_rhs_compatible_at_1e-8"]


def test_truncated_and_garbage_capture_rejected(tmp_path):
    path = fixture(tmp_path, [[2, -1], [-1, 2]], [1, 1])
    raw = path.read_bytes()
    path.write_bytes(raw[:-2])
    with pytest.raises(CertificateError, match="size"):
        analyze(path)
    path.write_bytes(b"WRONG" + raw[5:])
    with pytest.raises(CertificateError, match="magic"):
        analyze(path)


def test_dense_svd_bound_is_fail_closed(tmp_path):
    path = fixture(tmp_path, [[2, -1], [-1, 2]], [1, 1])
    with pytest.raises(CertificateError, match="dense SVD"):
        analyze(path, max_rows=1)


def test_ill_conditioned_but_invertible_not_null(tmp_path):
    path = fixture(tmp_path, [[1, 0], [0, 1e-8]], [1, 1])
    c = analyze(path)
    assert c["matrix"]["numerical_nullity_left"] == 0
    assert c["original_true_residual_relative_rhs"] < 1e-12
