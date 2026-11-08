"""Proof that identical total drag can hide different pressure and shear forces."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "scripts" / "naca_wall_radial_comparison.py"
SPEC = importlib.util.spec_from_file_location("naca_wall_radial_comparison", SOURCE)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def synthetic(pdrag: float, vdrag: float):
    rows = []
    for sign in (-1, 1):
        for i in range(64):
            x = (i + 0.5) / 64
            rows.append(
                {
                    "face": sign * 100 + i,
                    "x_over_c": x,
                    "y_over_c": sign * 0.06 * math.sin(math.pi * x),
                    "Cp": 1.0 + 0.1 * pdrag,
                    "Cf_tangent": 0.01 + 0.1 * vdrag,
                    "Cd_pressure_face": pdrag / 128,
                    "Cd_viscous_face": vdrag / 128,
                }
            )
    return rows


def cert(pdrag, vdrag):
    return {"cd": pdrag + vdrag, "cdp": pdrag, "cdv": vdrag}


def test_exact_total_cd_hides_large_opposing_component_errors():
    a = synthetic(0.04, 0.06)
    b = synthetic(0.03, 0.07)
    report = MODULE.compare(a, b, cert(0.04, 0.06), cert(0.03, 0.07))
    assert report["relative_total_cd_difference"] == pytest.approx(0.0)
    assert report["false_total_cd_agreement_risk"]
    assert report["relative_pressure_difference"] == pytest.approx(0.25)
    assert report["relative_viscous_difference"] == pytest.approx(1.0 / 7.0)
    assert report["force_cancellation_fraction"] == pytest.approx(1.0)
    assert report["gci"].startswith("NOT_PERFORMED")


def test_identical_profiles_do_not_trigger_cancellation():
    a = synthetic(0.04, 0.06)
    report = MODULE.compare(a, a, cert(0.04, 0.06), cert(0.04, 0.06))
    assert not report["false_total_cd_agreement_risk"]


def test_different_wall_geometry_fails_closed():
    a = synthetic(0.04, 0.06)
    b = synthetic(0.04, 0.06)
    b[0]["x_over_c"] += 0.001
    with pytest.raises(ValueError, match="geometrically comparable"):
        MODULE.compare(a, b, cert(0.04, 0.06), cert(0.04, 0.06))


def test_stale_force_certificate_fails_closed():
    a = synthetic(0.04, 0.06)
    with pytest.raises(ValueError, match="does not match face sum"):
        MODULE.compare(a, a, cert(0.05, 0.06), cert(0.04, 0.06))


def test_missing_face_fails_closed():
    a = synthetic(0.04, 0.06)
    with pytest.raises(ValueError, match="face counts"):
        MODULE.compare(a, a[:-1], cert(0.04, 0.06), cert(0.04, 0.06))
