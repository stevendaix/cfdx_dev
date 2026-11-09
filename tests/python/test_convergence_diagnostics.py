"""The monitoring math must never be misinterpreted as a convergence gate."""

import math

import pytest
from cfdx.convergence_diagnostics import (
    analyze_monitor_series,
    metric_group,
)
from cfdx.monitors import MonitorSample, MonitorSeries


def build_series(values, name="momentum_residual_relative"):
    return MonitorSeries(
        "solver",
        [
            MonitorSample(i, float(i), {name: value})
            for i, value in enumerate(values, 1)
        ],
    )


def test_descending_physical_residual_is_not_a_pass_certificate():
    result = analyze_monitor_series(build_series([1e-2, 1e-3, 1e-4, 1e-5, 1e-6]))
    assert result[0].group == "Nonlinear / equation"
    assert result[0].trend == "DECREASING"
    assert result[0].decades_per_iteration == pytest.approx(-1)
    assert result[0].source == "solver telemetry"
    assert not hasattr(result[0], "converged")


def test_stagnating_and_growing_and_missing_are_distinct():
    stagnant = analyze_monitor_series(build_series([1.0] * 6))[0]
    assert stagnant.trend == "STAGNATING"
    growing = analyze_monitor_series(build_series([0.001, 0.01, 0.1, 1]))[0]
    assert growing.trend == "GROWING"
    sparse = analyze_monitor_series(build_series([1, 0.1, 0.01]))[0]
    assert sparse.trend == "INSUFFICIENT"
    assert sparse.decades_per_iteration is None
    assert analyze_monitor_series(MonitorSeries("test", [])) == ()


def test_metric_categories_are_not_falsely_equated():
    assert metric_group("momentum_residual_relative") == "Nonlinear / equation"
    assert metric_group("continuity_normalized") == "Conservation"
    assert metric_group("flux_velocity_mismatch_linf") == "Conservation"
    assert metric_group("pressure_relative_residual") == "Linear solve"
    assert metric_group("pressure_linear_residual_true") == "Linear solve"
    assert metric_group("velocity_change_inf") == "State updates"
    assert metric_group("CFL") == "Other reported"


def test_missing_measurements_do_not_get_forward_filled():
    series = MonitorSeries(
        "solver",
        [
            MonitorSample(1, 0.0, {"momentum_residual_relative": 1.0}),
            MonitorSample(2, 0.0, {"continuity_normalized": 1e-3}),
            MonitorSample(3, 0.0, {"momentum_residual_relative": 0.1}),
            MonitorSample(4, 0.0, {"momentum_residual_relative": 0.01}),
        ],
    )
    result = {m.name: m for m in analyze_monitor_series(series)}
    assert result["momentum_residual_relative"].count == 3
    assert result["continuity_normalized"].count == 1
    assert all(row.trend == "INSUFFICIENT" for row in result.values())


def test_zero_and_signed_metrics_are_handled_without_log_nan():
    result = analyze_monitor_series(build_series([0, 0, -1e-20, 0]))[0]
    assert result.count == 4
    assert result.sparkline and len(result.sparkline) == 4
    assert math.isfinite(result.decades_per_iteration)


def test_non_finite_or_bad_limits_are_rejected():
    for bad in [float("nan"), float("inf"), True]:
        with pytest.raises(ValueError, match="non-finite"):
            analyze_monitor_series(build_series([bad]))
    with pytest.raises(ValueError, match="window"):
        analyze_monitor_series(build_series([]), window=2)
    with pytest.raises(ValueError, match="max_metrics"):
        analyze_monitor_series(build_series([]), max_metrics=0)


def test_latest_metric_slope_not_distorted_by_old_unrelated_cfl():
    series = MonitorSeries(
        "solver",
        [
            MonitorSample(
                i, float(i), {"CFL": 0.6, "momentum_residual_relative": 10 ** (-i)}
            )
            for i in range(1, 16)
        ],
    )
    rows = analyze_monitor_series(series, window=6)
    assert len(rows) == 1 and rows[0].count == 6
    assert rows[0].trend == "DECREASING"


def test_out_of_order_sample_is_not_silent_divergence():
    series = MonitorSeries(
        "solver",
        [
            MonitorSample(i, float(i), {"momentum_residual_relative": 10 ** (-i)})
            for i in (1, 3, 4, 8, 10)
        ],
    )
    assert analyze_monitor_series(series)[0].decades_per_iteration == pytest.approx(-1)
