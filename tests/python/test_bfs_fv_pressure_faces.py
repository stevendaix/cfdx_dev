"""Conservative implicit-FV-face bridge; never a full Rhie--Chow certificate."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from audit_bfs_fv_pressure_faces import (
    read_split,
    refine_implicit_faces,
    verify_face_operator,
)


def inputs():
    # Three cells, two internal oriented FV faces, two Dirichlet outlets.
    a = csr_matrix(np.array([[3, -2, 0], [-2, 4, -2], [0, -2, 3]], dtype=float))
    high = np.array([2.0, 3.0, 1.0])
    low = np.array([0.0, 1e-9, 0.0])
    rhs = np.asarray(a @ high + a @ low)
    owner = np.array([0, 1, 0, 2], dtype=np.int64)
    neighbour = np.array([1, 2, -1, -1], dtype=np.int64)
    outlet_faces = np.array([2, 3], dtype=np.int64)
    return a, rhs, high, low, owner, neighbour, outlet_faces


def test_real_fv_topology_conservative_split_operator():
    r = verify_face_operator(*inputs(), target_tolerance=1e-11)
    assert r["csr_face_topology_bijection"]
    assert r["real_mesh_implicit_face_gate_pass"]
    assert r["original_csr_linear_gate_pass"]
    assert abs(r["internal_face_pair_net"]) < 1e-18
    assert r["internal_fv_faces"] == 2
    assert r["outlet_pressure_faces"] == 2
    assert not r["physical_fv_face_flux_certified"]
    assert not r["actual_rhie_chow_predictor_and_nonorth_flux_replayed"]


def test_graph_mismatch_is_rejected():
    a, b, h, l, o, n, outlet = inputs()
    n[1] = 0
    with pytest.raises(
        ValueError, match="multiple physical faces|CSR off-diagonal graph"
    ):
        verify_face_operator(a, b, h, l, o, n, outlet, target_tolerance=1e-11)


def test_outlet_anchor_mismatch_is_rejected():
    a, b, h, l, o, n, outlet = inputs()
    with pytest.raises(ValueError, match="anchor"):
        verify_face_operator(a, b, h, l, o, n, outlet[:1], target_tolerance=1e-11)


def test_asymmetric_operator_is_rejected():
    a, b, h, l, o, n, outlet = inputs()
    a = a.tolil()
    a[1, 0] = -1.5
    with pytest.raises(ValueError, match="symmetric"):
        verify_face_operator(a.tocsr(), b, h, l, o, n, outlet, target_tolerance=1e-11)


def test_split_rejects_stale_checksum_and_invalid_rows(tmp_path):
    path = tmp_path / "p.split"
    path.write_text(
        "CFDX_PRESSURE_SPLIT_V1\n2 abcdef1234567890\n0x1p+0 0x0p+0\n0x2p+0 0x0p+0\n"
    )
    with pytest.raises(ValueError, match="different captured CSR"):
        read_split(path, 2, "deadbeef12345678")
    high, low = read_split(path, 2, "abcdef1234567890")
    assert np.array_equal(high, [1.0, 2.0])
    assert np.array_equal(low, [0.0, 0.0])
    path.write_text(path.read_text() + "0x1p+0 0x0p+0\n")
    with pytest.raises(ValueError, match="extra"):
        read_split(path, 2, "abcdef1234567890")


def test_zero_refinement_reports_both_gates_and_no_physics():
    r = refine_implicit_faces(*inputs(), target_tolerance=1e-11, iterations=0)
    assert r["final_implicit_face_gate_pass"]
    assert not r["physical_fv_flux_qualified"]
    assert r["no_production_state_or_rhie_chow_flux_modified"]
    with pytest.raises(ValueError, match="budget"):
        refine_implicit_faces(*inputs(), target_tolerance=1e-11, iterations=4)


def test_native_capture_checksum_is_recomputed_not_trusted(tmp_path):
    import struct

    from audit_bfs_fv_pressure_faces import verify_capture_checksum
    from seed_bfs_pressure_from_capture import read_capture

    data = [2.0, -1.0, -1.0, 2.0]
    columns = [0, 1, 0, 1]
    rhs = [1.0, 2.0]
    guess = [0.0, 0.0]
    tolerance = 1e-11
    bits = lambda x: struct.unpack("<Q", struct.pack("<d", x))[0]
    words = [2, 2, 4, 1, 0, 0, 100, bits(tolerance), 0, 2, 4]
    for col, value in zip(columns, data):
        words.extend([col, bits(value)])
    for b, x in zip(rhs, guess):
        words.extend([bits(b), bits(x)])
    raw = b"".join(int(word).to_bytes(8, "little") for word in words)
    state = 14695981039346656037
    for byte in raw:
        state = ((state ^ byte) * 1099511628211) & ((1 << 64) - 1)
    footer = f"{state:016x}"
    path = tmp_path / "case.cfdxlin"
    lines = ["CFDX_LINEAR_REPLAY_V1", "1 0 0", "2 2 4 100", f"{bits(tolerance):016x}"]
    lines.extend(f"{v:x}" for v in [0, 2, 4])
    lines.extend(f"{col:x} {bits(value):016x}" for col, value in zip(columns, data))
    lines.extend(f"{bits(b):016x} {bits(x):016x}" for b, x in zip(rhs, guess))
    lines.append(footer)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    a, b, _ = read_capture(path)
    verify_capture_checksum(path, a, b, footer)
    # A malicious/stale footer copied from the original must be rejected.
    path.write_text(
        path.read_text().replace(f"{bits(-1.0):016x}", f"{bits(-2.0):016x}", 1)
    )
    bad, bad_rhs, _ = read_capture(path)
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_capture_checksum(path, bad, bad_rhs, footer)
