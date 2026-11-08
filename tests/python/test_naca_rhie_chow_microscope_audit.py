from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT=Path(__file__).resolve().parents[2]/"scripts/naca_rhie_chow_microscope_audit.py"
spec=importlib.util.spec_from_file_location("naca_rhie_chow_microscope_audit",SCRIPT)
assert spec and spec.loader
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def sample(actual="0.012", expected="0.01", delta="0.002"):
    return ("=== CFDX MOMENTUM MICROSCOPE cell=201 cell_source=explicit algorithm=0 iteration=100 ===\n"
            "  face=5 owner=201 neighbour=200 patch=INTERNAL "
            f"rhie_chow_expected={expected} rhie_chow_actual={actual} delta_rhie_chow={delta}\n"
            "rhie_chow_expected_linf=0.01 rhie_chow_actual_linf=0.012 delta_rhie_chow_linf=0.002\n"
            "=== END CFDX MOMENTUM MICROSCOPE ===\n")


def test_parse_operator_vs_actual_is_diagnostic_only():
    r=module.parse_microscope(sample())
    assert r["cell"]==201 and r["iteration"]==100
    assert r["internal_face_count"]==1
    assert r["internal_max_absolute_difference"]==pytest.approx(0.002)
    assert "NOT_PHYSICAL" in r["status"]


@pytest.mark.parametrize("log",[
    "empty",
    sample().replace("=== END CFDX MOMENTUM MICROSCOPE ===",""),
    sample(delta="0.01"),
    sample(actual="nan"),
    sample().replace("neighbour=200","neighbour=-1"),
    sample().replace("delta_rhie_chow_linf=0.002","delta_rhie_chow_linf=0.001"),
])
def test_microscope_reader_fails_closed(log):
    with pytest.raises(ValueError):
        module.parse_microscope(log)
