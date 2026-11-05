"""V26 pressure-correction refinement must be explicit and fail closed."""
import json
from pathlib import Path

import h5py
import pytest

from scripts.generate_th03_openfoam_blockmesh import (
    make_case_setup, make_structured_mesh, write_case,
)


def test_v26_pressure_refinement_is_off_by_default():
    case = make_case_setup(1, 0.55, 0.3)
    assert "pressure_iterative_refinement_steps" not in case.numerics.raw_settings
    assert case.numerics.raw_settings["pressure_linear_tolerance"] == "1e-11"


def test_v26_pressure_refinement_exported_to_hdf5(tmp_path: Path):
    target = tmp_path / "refined.h5"
    write_case(target, make_structured_mesh("smoke"),
               mesh_source="v26-test", max_iterations=1,
               pressure_refinement_steps=4,
               pressure_tolerance=1e-10)
    with h5py.File(target) as f:
        values = json.loads(f.attrs["case_setup_json"])["numerics"]["raw_settings"]
    assert values["pressure_iterative_refinement_steps"] == "4"
    assert float(values["pressure_linear_tolerance"]) == 1e-10


@pytest.mark.parametrize("steps", [-1, 9, 99, 1.5])
def test_v26_reject_bad_pressure_refinement_step_counts(steps):
    with pytest.raises(ValueError, match="refinement steps"):
        make_case_setup(1, 0.55, 0.3, pressure_refinement_steps=steps)
