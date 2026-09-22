from cfdx.monitors import MonitorSample, MonitorSeries


def test_monitor_series_is_iteration_synchronized():
    series = MonitorSeries("pressure", [])
    series.append(MonitorSample(1, 0.1, {"p": 2.0}))
    series.append(MonitorSample(2, 0.2, {"p": 2.1}))
    assert series.at(2).time == 0.2


def test_monitor_series_upsert_merges_values_for_same_timestep():
    series = MonitorSeries("solver", [])
    series.upsert(MonitorSample(3, 0.3, {"CFL": 0.5}))
    series.upsert(MonitorSample(3, 0.3, {"p": 1.0e-3}))
    assert len(series.samples) == 1
    assert series.at(3).values == {"CFL": 0.5, "p": 1.0e-3}
