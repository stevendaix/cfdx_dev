"""Isolate geometric/mesh skew differences from airfoil contour changes."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))


def load(name: str):
    source = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GENERATOR = load("generate_naca_iso_aniso_gci")

arc_length_surface = GENERATOR.arc_length_surface
build_families = GENERATOR.build_families
load_generator = GENERATOR.load_generator
shear_interior_rings = GENERATOR.shear_interior_rings


@pytest.mark.parametrize("n", [64, 128])
def test_sheared_mesh_retains_exact_wall_and_farfield_nodes(n):
    generator = load_generator()
    original_surface = generator.naca0012_surface
    generator.naca0012_surface = lambda level, te=0.0: arc_length_surface(
        level, lambda dense: original_surface(dense, te)
    )
    nradial = 24 * n // 64
    nodes, cells, boundaries = generator.naca_o_grid(
        n, nradial, 20.0, radial_power=2.0, base_radial_layers=24
    )
    sheared_nodes, sheared_cells, sheared_boundary = shear_interior_rings(
        (nodes, cells, boundaries), n, nradial, strength=0.002
    )
    actual = np.asarray(sheared_nodes)
    original = np.asarray(nodes)
    circumferential = 2 * n + 1
    assert np.array_equal(actual[:circumferential], original[:circumferential])
    lo = nradial * circumferential
    hi = (nradial + 1) * circumferential
    assert np.array_equal(actual[lo:hi], original[lo:hi])
    assert not np.array_equal(actual[circumferential:lo], original[circumferential:lo])
    assert sheared_cells is cells
    assert sheared_boundary is boundaries


def test_shear_strength_is_physically_bounded(tmp_path):
    with pytest.raises(ValueError, match="strength"):
        build_families(
            tmp_path,
            levels=(64,),
            make_hdf5=False,
            only_family="sheared_isosurface",
            shear_strength=math.nan,
        )
    with pytest.raises(ValueError, match="strength"):
        build_families(
            tmp_path,
            levels=(64,),
            make_hdf5=False,
            only_family="sheared_isosurface",
            shear_strength=0.1,
        )


@pytest.mark.parametrize("family", ["quasi_isotropic", "anisotropic"])
def test_valid_n64_generation_end_to_end(family, tmp_path):
    rows = build_families(tmp_path, levels=(64,), make_hdf5=False, only_family=family)
    assert [row["n_surface"] for row in rows[family]] == [64]
    assert rows[family][0]["inverted_or_degenerate_cells"] == 0
    assert (tmp_path / family / "quality.json").exists()
