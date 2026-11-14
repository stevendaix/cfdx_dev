"""V18 wall owner inputs are explicit, fail closed, and keep the default intact."""

from __future__ import annotations

import pytest
from scripts.generate_th03_openfoam_blockmesh import make_case_setup


def test_v18_wall_owner_is_opt_in_only() -> None:
    base = make_case_setup(2, 0.25, 0.10)
    raw = base.numerics.raw_settings
    assert raw["v17_openfoam_wall_operator_experimental"] == "disabled"
    assert raw["v18_omega_wall_owner_constraints_experimental"] == "disabled"


def test_v18_requires_v17_atomic_wall_closure() -> None:
    with pytest.raises(ValueError, match="require V17"):
        make_case_setup(2, 0.25, 0.10, v18_omega_owner_cells=True)
    with pytest.raises(ValueError, match="kqR and molecular"):
        make_case_setup(2, 0.25, 0.10, v17_active_wall_operator=True)
    with pytest.raises(ValueError, match="kqR and molecular"):
        make_case_setup(2, 0.25, 0.10, v17_active_wall_operator=True,
                        openfoam_kqr_wall_v16=True)


def test_v18_stages_all_three_wall_conditions_together() -> None:
    setup = make_case_setup(
        2, 0.25, 0.10, v17_active_wall_operator=True,
        v18_omega_owner_cells=True, openfoam_kqr_wall_v16=True,
        sst_molecular_wall_diffusion=True,
    )
    raw = setup.numerics.raw_settings
    assert raw["v17_openfoam_wall_operator_experimental"] == "enabled"
    assert raw["v18_omega_wall_owner_constraints_experimental"] == "enabled"
    assert raw["sst_molecular_wall_diffusion"] == "true"
