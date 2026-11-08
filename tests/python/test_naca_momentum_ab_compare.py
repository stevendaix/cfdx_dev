from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "naca_momentum_ab_compare", SCRIPTS / "naca_momentum_ab_compare.py"
)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def histories() -> tuple[list[dict], list[dict]]:
    def one(value: float) -> list[dict]:
        return [{"iteration": i, "alpha_u": 0.6, "alpha_p": 0.25,
                 "high_order_fraction": 1.0,
                 **{key: value/(i+1) for key in module.METRICS}}
                for i in range(1, 11)]
    return one(3.0), one(1.0)


def budget(oscillating: bool) -> list[dict]:
    rows = []
    for i in range(1, 11):
        v = (1 if i % 2 == 0 else -1)*4.3e-5 if oscillating else 1e-4/i
        rows.append({
            "iteration": i, "component": 1, "rank": 0, "cell": 184,
            "x": -0.05, "y": 0.13, "z": 0.5,
            "convection": v + 1e-4, "viscous": 0.0,
            "transpose_viscous": 0.0, "pressure": -1e-4,
            "body_force": 0.0, "temporal": 0.0,
            "mass_divergence": 0.0, "bounded_divergence": 0.0,
            "collocated_difference": 1e-6,
            "physical_residual": v, "bounded_operator_residual": v,
            "worst_residual": abs(v), "algebraic_residual": v,
            "matrix_face_difference": 0.0, "matrix_face_scale": 0.5,
        })
    return rows


def test_matched_a_b_cycle_comparison_is_not_physical_certificate():
    tvd, up = histories()
    data = module.compare(tvd, up, budget(True), budget(False))
    assert data["accepted_iterations"] == 10
    assert data["tvd_transverse_hotspot"]["signed_term_cycle"]["signature"]["status"] == "PERIOD_TWO_SIGNATURE"
    assert data["final"]["momentum_rhs_relative"]["TVD_over_UPWIND"] == pytest.approx(3.0)
    assert data["tvd_convergence_proven"] is False


def test_rejects_unmatched_pressure_relaxation():
    tvd, up = histories()
    up[4]["alpha_p"] = 0.3
    with pytest.raises(ValueError, match="different"):
        module.compare(tvd, up, budget(True), budget(False))


def test_rejects_undocumented_continuation():
    tvd, up = histories()
    tvd[4]["high_order_fraction"] = 0.5
    with pytest.raises(ValueError, match="full high-order"):
        module.compare(tvd, up, budget(True), budget(False))


def test_rejects_budget_without_last_iteration():
    tvd, up = histories()
    with pytest.raises(ValueError, match="final accepted"):
        module.compare(tvd, up, budget(True)[:-1], budget(False))


@pytest.mark.parametrize("scheme,corruption", [
    ("TVD", "none"),
    ("UPWIND", "bad_scheme"),
    ("TVD", "nan"),
    ("TVD", "missing_iteration"),
])
def test_history_reader_fails_closed(tmp_path, scheme, corruption):
    path = tmp_path / "history.csv"
    headers = ["iteration", "convection_scheme", "alpha_u", "alpha_p",
               "high_order_fraction", *module.METRICS]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers)
        writer.writeheader()
        for i in range(1, 11):
            entry = dict.fromkeys(module.METRICS, 1e-4)
            entry.update(iteration=i, convection_scheme=scheme, alpha_u=0.6,
                         alpha_p=0.25, high_order_fraction=1.0)
            if corruption == "bad_scheme" and i == 5:
                entry["convection_scheme"] = "TVD"
            if corruption == "nan" and i == 6:
                entry["momentum_rhs_relative"] = "nan"
            if corruption == "missing_iteration" and i == 7:
                entry["iteration"] = "6"
            writer.writerow(entry)
    if corruption != "none":
        with pytest.raises(ValueError):
            module.read_history(path, scheme)
    else:
        assert len(module.read_history(path, scheme)) == 10
