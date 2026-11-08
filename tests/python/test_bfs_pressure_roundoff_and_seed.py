"""Independent proof tools for BFS pressure, not a solver qualification."""
from __future__ import annotations

import json
import struct
from pathlib import Path

import h5py
import numpy as np
import pytest

from scripts.audit_bfs_pressure_roundoff import audit, extended_residual
from scripts.seed_bfs_pressure_from_capture import build_seed, read_capture
from scripts.generate_th03_openfoam_blockmesh import make_structured_mesh, write_case


def _floathex(value: float) -> str:
    return f"{struct.unpack('=Q',struct.pack('=d', value))[0]:x}"


def capture(path: Path, n: int, bad=False) -> None:
    # A diagonal SPD matrix with fixed pressure reference and exact solution.
    with path.open("w") as s:
        s.write("CFDX_LINEAR_REPLAY_V1\n1 0 0\n")
        s.write(f"{n} {n} {n} 100\n{_floathex(1e-11)}\n")
        for i in range(n + 1):
            s.write(f"{i:x}\n")
        for i in range(n):
            s.write(f"{i:x} {_floathex(0.0 if bad and i == 0 else 2.0)}\n")
        for i in range(n):
            s.write(f"{_floathex(2.0 + i*0.01)} {_floathex(0.0)}\n")


def test_roundoff_audit_never_promotes_solver(tmp_path):
    p = tmp_path / "matrix.cfdxlin"
    capture(p, 5)
    result = audit(p, passes=2)
    assert result["status"] == "DIAGNOSTIC_ONLY_NOT_PHYSICAL_CERTIFICATE"
    assert result["rows"] == 5 and result["nonzeros"] == 5
    assert result["target_tolerance"] == 1e-11
    assert result["physical_convergence_proven"] is False
    assert result["steps"][-1]["extended_precision_true_relative"] < 1e-12


def test_bad_capture_and_singular_matrix_refused(tmp_path):
    cap = tmp_path / "a.cfdxlin"
    cap.write_text("INVALID\n")
    with pytest.raises(ValueError, match="not a CFDX pressure"):
        read_capture(cap)
    capture(cap, 5, bad=True)
    with pytest.raises(RuntimeError):
        audit(cap)


def test_seed_preserves_mesh_and_bc_and_is_not_qualified(tmp_path):
    case = tmp_path / "raw.cfdx.h5"
    output = tmp_path / "seeded.cfdx.h5"
    cap = tmp_path / "capture.cfdxlin"
    write_case(case, make_structured_mesh("smoke"), mesh_source="seed-test", max_iterations=1)
    with h5py.File(case) as f:
        n = int(f.attrs["n_cells"])
        original_setup = str(f.attrs["case_setup_json"])
        original_owner = f["owner"][:].copy()
    capture(cap, n)
    result = build_seed(case, cap, output)
    assert result["status"] == "DIAGNOSTIC_INITIAL_ONLY_NOT_QUALIFIED"
    assert result["pressure_rows"] == n
    with h5py.File(output) as f:
        assert f.attrs["case_setup_json"] == original_setup
        np.testing.assert_array_equal(f["owner"][:], original_owner)
        assert f["fields/scalar/p"].shape == (n,)
        assert f["fields/vector/U"].shape == (n, 3)
    with pytest.raises(ValueError, match="output would overwrite"):
        build_seed(case, cap, output)
    with pytest.raises(ValueError, match="blend"):
        build_seed(case, cap, tmp_path / "other.h5", blend=0.0)


def test_seed_rejects_incorrect_mesh_dimension(tmp_path):
    case = tmp_path / "raw.cfdx.h5"
    write_case(case, make_structured_mesh("smoke"), mesh_source="seed-test", max_iterations=1)
    cap = tmp_path / "capture.cfdxlin"
    capture(cap, 5)
    with pytest.raises(ValueError, match="does not match mesh"):
        build_seed(case, cap, tmp_path / "output.h5")
