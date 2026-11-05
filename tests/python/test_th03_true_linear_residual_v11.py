"""Independent pressure linear qualification: original matrix, not scaled system."""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

from scripts.audit_true_linear_residual_v11 import audit


def test_small_spd_direct_certification(tmp_path: Path):
    prefix = tmp_path / "a"
    m = csr_matrix(np.array([[2.0, -1.0], [-1.0, 2.0]]))
    mmwrite(str(prefix) + ".mtx", m)
    np.savetxt(str(prefix) + ".rhs", [1.0, 0.0])
    result = audit(prefix, target=1e-12)
    assert result["matrix_size"] == 2
    assert result["linear_status"] == "CERTIFIED"
    assert result["physical_status"] == "NOT_EVALUATED"
    assert result["normwise_backward_error"] < 1e-14


def test_invalid_input_is_rejected(tmp_path: Path):
    prefix = tmp_path / "a"
    mmwrite(str(prefix) + ".mtx", csr_matrix(np.eye(2)))
    np.savetxt(str(prefix) + ".rhs", [1.0, 1.0])
    with pytest.raises(ValueError):
        audit(prefix, target=0)
    with pytest.raises(ValueError):
        audit(prefix, target=float("nan"))


def test_nonsquare_fails_closed(tmp_path: Path):
    prefix = tmp_path / "a"
    mmwrite(str(prefix) + ".mtx", csr_matrix(np.ones((2, 3))))
    np.savetxt(str(prefix) + ".rhs", [1.0, 1.0])
    with pytest.raises(ValueError, match="invalid shape"):
        audit(prefix)
