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


def test_monitor_series_json_roundtrip_preserves_iteration_time_and_values(tmp_path):
    series = MonitorSeries("solver", [])
    series.upsert(MonitorSample(1, 0.1, {"CFL": 0.5}))
    series.upsert(MonitorSample(2, 0.2, {"p": 2.0e-3}))

    path = series.write_json(tmp_path / "monitor.json")
    restored = MonitorSeries.read_json(path)

    assert restored.name == "solver"
    assert [s.iteration for s in restored.samples] == [1, 2]
    assert [s.time for s in restored.samples] == [0.1, 0.2]
    assert restored.at(1).values["CFL"] == 0.5
    assert restored.at(2).values["p"] == 2.0e-3


def test_monitor_series_rejects_non_finite_persisted_values():
    payload = {
        "name": "solver",
        "samples": [{"iteration": 1, "time": 0.1, "values": {"p": float("nan")}}],
    }
    import pytest

    with pytest.raises(ValueError, match="invalid monitor value"):
        MonitorSeries.from_dict(payload)


def test_monitor_series_rejects_non_monotonic_persisted_samples():
    payload = {
        "name": "solver",
        "samples": [
            {"iteration": 2, "time": 0.2, "values": {"p": 1.0}},
            {"iteration": 1, "time": 0.1, "values": {"p": 2.0}},
        ],
    }
    import pytest

    with pytest.raises(ValueError, match="monotonic"):
        MonitorSeries.from_dict(payload)
