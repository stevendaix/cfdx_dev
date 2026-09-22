from pathlib import Path

import pytest

from cfdx.post import read_monitor_csv, reduce_values


def test_read_monitor_csv(tmp_path: Path) -> None:
    path = tmp_path / "forces.csv"
    path.write_text("time,drag,lift\n0,1,2\n1,3,4\n", encoding="utf-8")
    points = read_monitor_csv(path)
    assert points[1].time == 1
    assert points[1].values["drag"] == 3


def test_read_monitor_csv_requires_time(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text("iteration,value\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="time"):
        read_monitor_csv(path)


def test_reduce_values() -> None:
    assert reduce_values([1, 2, 5]) == {
        "count": 3.0,
        "min": 1.0,
        "max": 5.0,
        "mean": 8 / 3,
    }
