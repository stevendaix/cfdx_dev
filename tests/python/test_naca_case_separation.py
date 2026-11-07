"""Regressions: a TMR grid is not a SST aerodynamic validation."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "naca_case_separation", ROOT / "scripts" / "validate_naca_case_separation.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def registry():
    return json.loads(
        (ROOT / "docs" / "validation" / "CFDX_QUALIFICATION_REGISTRY.json").read_text(
            encoding="utf-8"
        )
    )


def test_separate_existing_laminar_and_turbulent_cases():
    laminar, turbulent = MODULE.check_registry(registry())
    assert laminar["id"] == "EXT-NACA0012"
    assert turbulent["id"] == "TUR-NACA0012"
    assert turbulent["status"] != "PASS"


@pytest.mark.parametrize(
    "mutator",
    [
        lambda s: s.update(status="PASS"),
        lambda s: s.update(numerical_schemes="laminar"),
        lambda s: s.update(parameters="Re=1000; alpha=0"),
    ],
)
def test_sst_cannot_be_relabelled_or_promoted_by_metadata(mutator):
    data = copy.deepcopy(registry())
    sst = next(c for c in data["cases"] if c["id"] == "TUR-NACA0012")
    mutator(sst)
    with pytest.raises(ValueError):
        MODULE.check_registry(data)


def test_unproven_tmr_mesh_is_not_automatically_sst(tmp_path):
    import h5py
    import numpy as np

    path = tmp_path / "tmr_mesh.h5"
    with h5py.File(path, "w") as h:
        h.attrs["validation_source"] = "NASA TMR NACA0012"
        h.attrs["validation_status"] = "imported-provenance-only"
        h.attrs["validation_source_url"] = (
            "https://tmbwg.github.io/turbmodels/naca0012_grids.html"
        )
        h.attrs["validation_source_sha256"] = "a" * 64
        h.attrs["validation_grid"] = "225x65"
        h.attrs["boundary_patches"] = "airfoil:0:1:0;front:1:1:6;back:2:1:6"
        h.create_dataset("points", data=np.zeros((8, 3)))
        h.create_dataset("owner", data=np.zeros(3, dtype=np.uint64))
        h.create_dataset("neighbour", data=-np.ones(3, dtype=np.int64))
        h.create_dataset("cell_offsets", data=np.array([0, 6], dtype=np.uint64))
    assert MODULE.check_tmr_mesh(path) == "225x65"
    with h5py.File(path, "r+") as h:
        h.attrs["validation_status"] = "PASS"
    with pytest.raises(ValueError, match="import-only"):
        MODULE.check_tmr_mesh(path)


def test_reject_generated_laminar_mesh_as_tmr(tmp_path):
    import h5py

    path = tmp_path / "laminar_ogrid.h5"
    with h5py.File(path, "w") as h:
        h.attrs["validation_source"] = "CFDX generated O-grid"
    with pytest.raises(ValueError, match="NASA/TMR"):
        MODULE.check_tmr_mesh(path)
