"""Pressure audit remains a separate independent numerical certificate."""
from pathlib import Path

import numpy as np
import pytest
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from scripts.audit_pressure_precision_v15 import inspect


def test_well_conditioned_pressure_matrix_is_certified(tmp_path: Path):
    matrix = tmp_path / 'A.mtx'
    rhs = tmp_path / 'b.txt'
    mmwrite(matrix, csr_matrix(np.array([[4.,1.], [1.,3.]])))
    np.savetxt(rhs, [1.,2.])
    result = inspect(matrix, rhs, 1e-11, 2)
    assert result['n'] == 2
    assert result['best_true_residual_fp64'] < 1e-11
    assert result['physical_qualification'] == 'NOT_QUALIFIED'


def test_pressure_audit_rejects_invalid_tolerance(tmp_path: Path):
    with pytest.raises(ValueError, match='invalid pressure audit'):
        inspect(tmp_path / 'none.mtx', tmp_path / 'none.rhs', 0.)