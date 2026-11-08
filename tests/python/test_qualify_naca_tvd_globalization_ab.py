import importlib.util
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
FILE = SCRIPTS / "qualify_naca_tvd_globalization_ab.py"
spec = importlib.util.spec_from_file_location("picard_ab", FILE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def line(m=1e-9, active=1e-9, cont=1e-9):
    return (
        "INCOMPRESSIBLE_ITER iter=5 momentum_residual_rhs_relative="
        + str(m)
        + " momentum_active_component_residual_max="
        + str(active)
        + " continuity_normalized="
        + str(cont)
        + " pressure_relative_residual=1e-12 velocity_change_inf=1e-10"
        + " pressure_change_inf=1e-10 high_order_fraction=1 full_order_eligible=1"
    )


def test_good_only_iterative_not_reference():
    x = module.diagnose("NACA0012_SINGLE: PASS", line(), 0, "baseline")
    assert x["solver_contract_pass"]
    assert not x["physical_reference_qualified"]
    assert x["reason"] == "ITERATIVE_CONVERGENCE_ONLY"


@pytest.mark.parametrize(
    "stderr,stdout,exit_code",
    [
        (line(m=1e-4), "NACA0012_SINGLE: PASS", 0),
        (line(active=1e-4), "NACA0012_SINGLE: PASS", 0),
        (line(cont=1e-5), "NACA0012_SINGLE: PASS", 0),
        (line(), "NACA0012_SINGLE: PASS", 1),
        (line(), "", 0),
        (
            line().replace("high_order_fraction=1", "high_order_fraction=0.8"),
            "NACA0012_SINGLE: PASS",
            0,
        ),
        (
            line() + "\nNACA0012_QUALIFICATION: FAIL: fake pass",
            "NACA0012_SINGLE: PASS",
            0,
        ),
        (
            "INCOMPRESSIBLE_ITER iter=1 momentum_residual_rhs_relative=1e-12",
            "NACA0012_SINGLE: PASS",
            0,
        ),
    ],
)
def test_fail_closed(stderr, stdout, exit_code):
    verdict = module.diagnose(stdout, stderr, exit_code, "baseline")
    assert not verdict["solver_contract_pass"]


def test_memory_requires_source_metric():
    a = module.diagnose("NACA0012_SINGLE: PASS", line(), 0, "memory_raw_025")
    assert not a["solver_contract_pass"]
    b = module.diagnose(
        "NACA0012_SINGLE: PASS",
        line() + "\nTVD_DEFECT_MEMORY source_relative=3e-5",
        0,
        "memory_raw_025",
    )
    assert not b["solver_contract_pass"]
    c = module.diagnose(
        "NACA0012_SINGLE: PASS",
        line() + "\nTVD_DEFECT_MEMORY source_relative=1e-10",
        0,
        "memory_raw_025",
    )
    assert c["solver_contract_pass"]
