from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def load(name: str):
    source = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CROSS = load("naca_cross_surface_profile_comparison")


def test_face_bin_deltas_conserve_components():
    a = {
        "force": {"cd": 0.11, "cdp": 0.04, "cdv": 0.07, "face_count": 2},
        "bins": [
            {
                "side": "upper",
                "range": [0.0, 1.0],
                "cdp": 0.02,
                "cdv": 0.04,
                "face_count": 2,
            },
            {
                "side": "lower",
                "range": [0.0, 1.0],
                "cdp": 0.02,
                "cdv": 0.03,
                "face_count": 2,
            },
        ],
    }
    b = {
        "force": {"cd": 0.12, "cdp": 0.038, "cdv": 0.082},
        "bins": [
            {
                "side": "upper",
                "range": [0.0, 1.0],
                "cdp": 0.019,
                "cdv": 0.046,
                "face_count": 2,
            },
            {
                "side": "lower",
                "range": [0.0, 1.0],
                "cdp": 0.019,
                "cdv": 0.036,
                "face_count": 2,
            },
        ],
    }
    r = CROSS.compare(a, b)
    assert r["delta_force"]["cdp"] == pytest.approx(-0.002)
    assert r["delta_force"]["cdv"] == pytest.approx(0.012)
    assert r["delta_force"]["cd"] == pytest.approx(0.01)
    assert r["pressure_viscous_cancellation_fraction"] == pytest.approx(
        1 - 0.01 / 0.014
    )


def test_nonconservative_bin_rejected():
    a = {
        "force": {"cd": 0.11, "cdp": 0.04, "cdv": 0.07, "face_count": 2},
        "bins": [
            {
                "side": "upper",
                "range": [0.0, 1.0],
                "cdp": 0.04,
                "cdv": 0.07,
                "face_count": 2,
            }
        ],
    }
    b = {
        "force": {"cd": 0.12, "cdp": 0.038, "cdv": 0.082},
        "bins": [
            {
                "side": "upper",
                "range": [0.0, 1.0],
                "cdp": 0.040,
                "cdv": 0.082,
                "face_count": 2,
            }
        ],
    }
    with pytest.raises(ValueError, match="integrated bin"):
        CROSS.compare(a, b)
