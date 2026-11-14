"""True CSR linear qualification must not be mistaken for physical validation."""
from pathlib import Path
import numpy as np
import pytest
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from scripts.qualify_coupled_matrix_v16 import qualify


def test_independent_both_blocks_converge_without_nonphysical_pass(tmp_path:Path):
    prefix=tmp_path/'matrix'
    A=csr_matrix([[2.,0.,0.,1.],[0.,3.,0.,1.],
                  [0.,0.,4.,-1.],[1.,1.,-1.,0.]])
    mmwrite(str(prefix)+'_A.mtx',A)
    mmwrite(str(prefix)+'_b.mtx',np.array([[1.],[2.],[3.],[1.]]))
    mmwrite(str(prefix)+'_x0.mtx',np.zeros((4,1)))
    result=qualify(prefix)
    assert result['direct_lu']['full_relative']<1e-12
    assert result['ilut_gmres']['original_system_converged']
    assert result['physical_validation']=='NOT_QUALIFIED'


def test_missing_rhs_is_fail_closed(tmp_path:Path):
    with pytest.raises(ValueError,match='capture'):
        qualify(tmp_path/'nothing')
