"""Protect the existing V4 OpenFOAM six-block BFS, not a replacement mesh.

This test is geometry/provenance ONLY. OpenFOAM polyMesh identity and
Driver-Seegmiller physical qualification require separate evidence.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import generate_th03_openfoam_blockmesh as reference


def test_openfoam_v2606_six_block_native_geometry_is_preserved():
    specs = reference.block_specs("full")
    assert len(specs) == 6
    assert [(x["nx"], x["ny"]) for x in specs] == [
        (2, 65),
        (88, 65),
        (97, 48),
        (33, 48),
        (97, 65),
        (33, 65),
    ]
    assert sum(x["nx"] * x["ny"] for x in specs) == 20540
    assert specs[0]["x0"] == -130 and specs[-1]["x1"] == 50
    assert reference.H == pytest.approx(0.0127)
    assert tuple(reference.PATCH_TYPES) == (
        "inlet",
        "outlet",
        "lowerWallStartup",
        "upperWallStartup",
        "upperWall",
        "lowerWall",
        "front",
        "back",
    )
    assert reference.PATCH_TYPES["front"] == "empty"
    assert reference.PATCH_TYPES["lowerWallStartup"] == "symmetryPlane"


def test_openfoam_tutorial_and_driver_seegmiller_are_distinct_contracts():
    # CFDX's V4 generator targets the OpenFOAM tutorial's nominal Re=36k.
    # TH-03 frozen experimental spec sets 37,500, a DIFFERENT viscosity.
    reynolds_tutorial = reference.VELOCITY * reference.H / reference.NU
    assert reynolds_tutorial == pytest.approx(35983.333333333336, rel=1e-6)
    nu_driver = reference.VELOCITY * reference.H / 37500.0
    assert nu_driver != pytest.approx(reference.NU, rel=1e-3)
    assert reference.OMEGA_IN == pytest.approx(
        reference.K_IN / (reference.NU * reference.INLET_TURBULENT_VISCOSITY_RATIO)
    )
