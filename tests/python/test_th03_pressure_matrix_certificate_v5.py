"""Pressure matrix independent certification fails closed on the original residual."""

from pathlib import Path

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

from scripts.analyse_th03_pressure_v5 import analyse


def test_exact_small_pressure_poisson_matrix(tmp_path: Path):
    p = tmp_path / "small"
    m = csr_matrix(np.array([[2.0, -1.0], [-1.0, 2.0]]))
    mmwrite(str(p) + ".mtx", m)
    np.savetxt(str(p) + ".rhs", [1.0, 0.0])
    report = analyse(p)
    assert report["rows"] == 2
    assert report["symmetric_max_absolute_error"] == 0
    assert report["status"] == "LINEAR_CERTIFIED"


def test_stricter_than_fp64_roundoff_fails_closed(tmp_path: Path):
    p = tmp_path / "strict"
    rng = np.random.default_rng(7)
    m = rng.standard_normal((7, 7))
    mmwrite(str(p) + ".mtx", csr_matrix(m.T @ m + 0.1 * np.eye(7)))
    np.savetxt(str(p) + ".rhs", rng.standard_normal(7))
    report = analyse(p, 1e-35)
    assert report["status"] == "LINEAR_NOT_CERTIFIED"
    assert report["direct_solver_true_relative_residual"] > 1e-35
