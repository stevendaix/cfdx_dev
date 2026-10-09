"""Regression tests for ambiguous probe CSV histories and metadata."""

import pytest
from cfdx.probe_validation import (
    PROBE_CSV_MAGIC,
    PROBE_CSV_VERSION,
    ProbeSample,
    ProbeSeries,
    read_probe_csv,
    write_probe_csv,
)


def test_write_probe_csv_rejects_partial_metadata(tmp_path):
    series = ProbeSeries(
        "pressure",
        (0.5, 0.5, 0.5),
        (ProbeSample(1, 0.1, 1.0),),
        field="p",
    )
    with pytest.raises(ValueError, match="incomplete metadata"):
        write_probe_csv(tmp_path / "partial_metadata.csv", [series])


def test_write_probe_csv_rejects_duplicate_iterations(tmp_path):
    series = ProbeSeries(
        "pressure",
        (),
        (ProbeSample(1, 0.1, 1.0), ProbeSample(1, 0.2, 1.1)),
    )
    with pytest.raises(ValueError, match="duplicate iteration 1"):
        write_probe_csv(tmp_path / "duplicate_iteration.csv", [series])


def test_read_probe_csv_rejects_duplicate_iterations(tmp_path):
    path = tmp_path / "duplicate_iteration.csv"
    path.write_text(
        f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n"
        "# columns: probe,iteration,time,value\n"
        "pressure,1,0.1,1.0\n"
        "pressure,1,0.2,1.1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="strictly increasing iterations"):
        read_probe_csv(path)


def test_write_probe_csv_rejects_duplicate_probe_names(tmp_path):
    series = [
        ProbeSeries("pressure", (), (ProbeSample(1, 0.1, 1.0),)),
        ProbeSeries("pressure", (), (ProbeSample(2, 0.2, 1.1),)),
    ]
    with pytest.raises(ValueError, match="probe names must be unique"):
        write_probe_csv(tmp_path / "duplicate_probe.csv", series)


def test_write_probe_csv_creates_parent_directory(tmp_path):
    destination = tmp_path / "results" / "probes" / "monitor.csv"
    series = [ProbeSeries("pressure", (), (ProbeSample(1, 0.1, 1.0),))]
    assert write_probe_csv(destination, series) == destination
    assert destination.is_file()



def test_write_probe_csv_rejects_decreasing_time_in_iteration_order(tmp_path):
    series = ProbeSeries(
        "pressure",
        (),
        (ProbeSample(1, 0.2, 10.0), ProbeSample(2, 0.1, 11.0)),
    )
    with pytest.raises(ValueError, match="non-decreasing time"):
        write_probe_csv(tmp_path / "decreasing_time.csv", [series])


def test_write_probe_csv_accepts_constant_time(tmp_path):
    destination = tmp_path / "constant_time.csv"
    series = ProbeSeries(
        "pressure",
        (),
        (ProbeSample(1, 0.1, 10.0), ProbeSample(2, 0.1, 11.0)),
    )
    write_probe_csv(destination, [series])
    assert [sample.time for sample in read_probe_csv(destination)[0].samples] == [
        0.1,
        0.1,
    ]


def test_write_probe_csv_allows_shared_iterations_across_probes(tmp_path):
    destination = tmp_path / "shared_iterations.csv"
    series = [
        ProbeSeries("pressure", (), (ProbeSample(1, 0.1, 10.0),)),
        ProbeSeries("velocity", (), (ProbeSample(1, 0.1, 2.0),)),
    ]
    write_probe_csv(destination, series)
    restored = {probe.name: probe for probe in read_probe_csv(destination)}
    assert restored["pressure"].samples[0].iteration == 1
    assert restored["velocity"].samples[0].iteration == 1
