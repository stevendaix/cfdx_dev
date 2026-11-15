"""No algebraic compensated result can be promoted to physical BFS success."""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_bfs_compensated_pressure import compensated_residual, inspect, two_sum_update


def test_two_sum_survives_fp64_cancellation():
    hi = np.array([1e8, 1e8], dtype=np.float64)
    lo = np.zeros(2, dtype=np.float64)
    two_sum_update(hi, lo, np.array([0.0, 2e-9]))
    assert hi[1] == 1e8
    assert lo[1] == 2e-9
    a = csr_matrix([[2.0, -2.0], [-2.0, 2.0]])
    b = np.array([-4e-9, 4e-9])
    residual, _, graph = compensated_residual(a, b, hi, lo)
    assert graph
    assert float(np.linalg.norm(residual)) < 1e-20
    assert float(np.linalg.norm(b - a @ hi)) > 1e-9


def test_spd_anchored_system_has_true_extended_residual_and_no_fv_promotion():
    a = csr_matrix([[4.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 3.0]])
    b = np.array([1.0, 2.0, 3.0])
    report = inspect(a, b, tolerance=1e-12, refinements=1)
    assert report["observations"][-1]["compensated_algebraic_gate_passed"]
    assert report["pressure_operator_symmetric_m_matrix"]
    assert report["algebraic_face_audit"]["algebraic_face_relative_residual"] < 1e-12
    assert (
        abs(report["algebraic_face_audit"]["internal_face_flux_pair_global_sum"])
        < 1e-15
    )
    assert not report["physical_face_flux_certified"]
    assert not report["accept_as_production_pressure_correction"]


def test_nonsymmetric_matrix_is_not_graph_certified():
    a = csr_matrix([[2.0, -1.0], [-2.0, 3.0]])
    report = inspect(a, np.array([1.0, 2.0]), refinements=1)
    assert not report["pressure_operator_symmetric_m_matrix"]
    assert report["algebraic_face_audit"] is None
    assert not report["physical_face_flux_certified"]


@pytest.mark.parametrize("invalid", [0.0, float("nan"), -1.0])
def test_bad_precision_gate_is_rejected(invalid):
    with pytest.raises(ValueError):
        inspect(csr_matrix([[2.0]]), np.array([1.0]), tolerance=invalid)


def test_homogeneous_rhs_is_not_a_false_success():
    with pytest.raises(ValueError):
        inspect(csr_matrix([[2.0]]), np.array([0.0]))


def test_actual_bfs_capture_if_explicitly_supplied():
    import os

    cap = os.environ.get("CFDX_BFS_PRESSURE_CAPTURE")
    if not cap:
        pytest.skip("requires immutable captured BFS V4 pressure matrix")
    from seed_bfs_pressure_from_capture import read_capture

    a, b, metadata = read_capture(Path(cap))
    report = inspect(a, b, tolerance=metadata["target_tolerance"], refinements=1)
    assert report["n"] == 20540
    assert report["nnz"] == 102034
    assert report["pressure_operator_symmetric_m_matrix"]
    assert not report["physical_face_flux_certified"]
