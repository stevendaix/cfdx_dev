from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from analyse_thermohydraulic_timeseries import TimeseriesError, compute


def make_series(path: Path, n=1000, thermal=True, drift=False, missing=False):
    fields = ["time_s", "station", "u_x_m_s", "u_y_m_s", "u_z_m_s"]
    if thermal:
        fields += ["temperature_k"]
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        for j in range(n):
            time = j * 0.02
            z = math.sin(2 * math.pi * 5 * time)
            for station in ("wall", "core"):
                if missing and j == 50 and station == "core":
                    continue
                row = {
                    "time_s": time,
                    "station": station,
                    "u_x_m_s": 2 + z + (j / n if drift else 0),
                    "u_y_m_s": z * 0.5,
                    "u_z_m_s": z * 0.2,
                }
                if thermal:
                    row["temperature_k"] = 310 + z
                writer.writerow(row)


def test_known_sinusoidal_covariance_heat_flux_and_frequency(tmp_path):
    path = tmp_path / "sine.csv"
    make_series(path)
    result = compute(path, min_effective=2)
    s = result["stations"]["core"]
    assert s["mean"]["u_x_m_s"] == pytest.approx(2, abs=1e-12)
    assert s["reynolds_uv_m2_s2"] == pytest.approx(0.25, rel=0.003)
    assert s["turbulent_heat_flux_uT_mK_s"] == pytest.approx(0.5, rel=0.003)
    f = np.array(s["spectrum_u_x"]["frequencies_hz"])
    power = np.array(s["spectrum_u_x"]["psd_m2_s"])
    assert f[np.argmax(power)] == pytest.approx(5, abs=0.05)
    assert result["physical_qualification"] == "NOT_QUALIFIED"


def test_nonuniform_and_missing_probe_rejected(tmp_path):
    path = tmp_path / "missing.csv"
    make_series(path, missing=True)
    with pytest.raises(TimeseriesError, match="exactly the same"):
        compute(path)
    path = tmp_path / "short.csv"
    make_series(path, n=9)
    with pytest.raises(TimeseriesError, match="not enough"):
        compute(path)


def test_mean_drift_detected(tmp_path):
    path = tmp_path / "drift.csv"
    make_series(path, drift=True)
    data = compute(path, min_effective=2, max_relative_drift=0.02)
    assert data["status"] == "INCOMPLETE_TIME_STATISTICS"


def test_no_temperature_is_valid_statistics_not_thermal_validation(tmp_path):
    path = tmp_path / "velocity.csv"
    make_series(path, thermal=False)
    data = compute(path, min_effective=2)
    assert "turbulent_heat_flux_uT_mK_s" not in data["stations"]["core"]
    assert data["physical_qualification"] == "NOT_QUALIFIED"
