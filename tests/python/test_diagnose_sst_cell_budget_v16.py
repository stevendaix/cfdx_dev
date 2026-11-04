"""V16 row-budget parser rejects incomplete rows and accepts signed conservation."""
from pathlib import Path
import pytest
from scripts.diagnose_sst_cell_budget_v16 import inspect


def test_budget_parser_round_trip(tmp_path: Path):
    log = tmp_path / 'budget.log'
    log.write_text(
        'SST_CELL_BUDGET iteration=2 field=k cell=210 convection=1.0 '
        'diffusion=2.0 production=-1.5 destruction=1.0 sum=2.5 residual=2.5 '
        'decomposition_error=0 source_identity_error=0 physical_production=-1.5 '
        'physical_destruction=1 physical_cross=0 local_relative=0.1 '
        'div_phi=0.1 wall_faces=1 wall_orth_flux=2.5 '
        'wall_owner_flux=5 wall_molecular_flux=2.5 qualified=0\n'
        'SST_CERTIFIED_INNER iteration=2\nConverged NO Iterations 3\n')
    result=inspect(log)
    assert result['sample_count']==1
    assert result['samples'][0]['wall_owner_molecular_ratio']==2.0
    assert result['terminal']=={'converged':False,'outer_iterations':3}
    assert result['internal_sst_certifications']==1


def test_budget_parser_rejects_wrong_sum_and_nonfinite(tmp_path: Path):
    log=tmp_path/'bad.log'
    line=('SST_CELL_BUDGET iteration=1 field=k cell=1 convection=1 '
          'diffusion=0 production=0 destruction=0 sum=1 residual=1 '
          'decomposition_error=0 source_identity_error=0 physical_production=0 '
          'physical_destruction=0 physical_cross=0 local_relative=0.1 '
          'div_phi=0 wall_faces=0 wall_orth_flux=0 '
          'wall_owner_flux=0 wall_molecular_flux=0 qualified=0\n')
    for corrupt in (line.replace('sum=1', 'sum=2'),line.replace('residual=1','residual=nan'),
                    line.replace('qualified=0','qualified=1')):
        log.write_text(corrupt)
        with pytest.raises(ValueError):
            inspect(log)