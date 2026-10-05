import math
import textwrap

import pytest

from cfdx.probe_validation import (
    PROBE_CSV_MAGIC,
    PROBE_CSV_VERSION,
    ProbeSample,
    ProbeSeries,
    read_probe_csv,
    write_probe_csv,
)

def test_probe_series_validates_quantitatively():
    series=ProbeSeries("pressure",(0.5,0.5,0.5),(
        ProbeSample(1,0.1,1.0),ProbeSample(2,0.2,1.1)))
    assert series.validate_against((1.0,1.1),rtol=1e-12)==0.0

def test_probe_series_rejects_out_of_tolerance():
    series=ProbeSeries("pressure",(0.5,), (ProbeSample(1,0.1,1.2),))
    with pytest.raises(ValueError,match="exceeds tolerance"):
        series.validate_against((1.0,),rtol=1e-8)

def test_probe_series_rejects_non_monotonic_samples():
    series=ProbeSeries("pressure",(0.5,),(
        ProbeSample(2,0.2,1.0),ProbeSample(1,0.1,1.0)))
    with pytest.raises(ValueError,match="monotonic"):
        series.validate_against((1.0,1.0))


def test_write_probe_csv_emits_native_schema(tmp_path):
    path = write_probe_csv(tmp_path / "probes.csv", [
        ProbeSeries("pressure", (), (ProbeSample(1, 0.1, 1.0),)),
    ])
    text = path.read_text()
    assert text.startswith(f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n")
    assert text.splitlines()[1] == "# columns: probe,iteration,time,value"
    assert text.splitlines()[2] == "pressure,1,0.1,1.0"


def test_write_probe_csv_byte_stable(tmp_path):
    # repr() yields the shortest float that round-trips, the same spelling the
    # C++ writer emits via std::to_chars, so the export is byte-stable across the
    # two implementations for any value with a non-integer magnitude.
    path = write_probe_csv(tmp_path / "probes.csv", [
        ProbeSeries("u_mag", (), (
            ProbeSample(1, 1.0, 0.0),
            ProbeSample(2, 1.0 / 3.0, 1.0 / 3.0),
            ProbeSample(3, 1e-9, 123.45678901234568),
        )),
    ])
    expected = (
        f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n"
        "# columns: probe,iteration,time,value\n"
        f"u_mag,1,{repr(1.0)},{repr(0.0)}\n"
        f"u_mag,2,{repr(1.0 / 3.0)},{repr(1.0 / 3.0)}\n"
        f"u_mag,3,{repr(1e-9)},{repr(123.45678901234568)}\n"
    )
    assert path.read_text() == expected


def test_write_probe_csv_orders_by_name_then_iteration(tmp_path):
    path = write_probe_csv(tmp_path / "probes.csv", [
        ProbeSeries("zeta", (), (ProbeSample(2, 0.2, 3.0),)),
        ProbeSeries("alpha", (), (
            ProbeSample(2, 0.2, 2.0),
            ProbeSample(1, 0.1, 1.0),
        )),
    ])
    rows = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    assert rows == [
        "alpha,1,0.1,1.0",
        "alpha,2,0.2,2.0",
        "zeta,2,0.2,3.0",
    ]


def test_read_probe_csv_round_trips(tmp_path):
    series = [
        ProbeSeries("pressure", (), (
            ProbeSample(1, 0.1, 1.0),
            ProbeSample(2, 0.2, 1.1),
        )),
        ProbeSeries("u_x", (), (ProbeSample(1, 0.1, 0.5),)),
    ]
    path = write_probe_csv(tmp_path / "probes.csv", series)
    loaded = read_probe_csv(path)
    assert tuple(s.name for s in loaded) == ("pressure", "u_x")
    for original, read_back in zip(series, loaded):
        assert original.name == read_back.name
        assert read_back.point == ()
        assert read_back.samples == original.samples


def test_read_probe_csv_parses_native_c_format(tmp_path):
    # Hand-written bytes in the exact format cfdx::io::write_probe_csv produces:
    # shortest round-trip floats with a decimal point, max_iterations ordering.
    csv = tmp_path / "from_cpp.csv"
    csv.write_text(textwrap.dedent("""\
        # cfdx-probe-csv v1
        # columns: probe,iteration,time,value
        pressure,1,0.5,0.3333333333333333
        pressure,2,0.5,0.3333333333333333
        u_mag,1,0.5,0.0
    """))
    loaded = read_probe_csv(csv)
    assert len(loaded) == 2
    pressure, u_mag = loaded
    assert pressure.name == "pressure"
    assert len(pressure.samples) == 2
    assert pressure.samples[0].value == pytest.approx(1.0 / 3.0)
    assert u_mag.samples[0].value == 0.0


@pytest.mark.parametrize("bad", [
    "pressure",
    "pressure,1,0,0.1,x",
])
def test_read_probe_csv_rejects_bad_row(tmp_path, bad):
    csv = tmp_path / "bad.csv"
    csv.write_text(f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n# columns: probe,iteration,time,value\n{bad}\n")
    with pytest.raises(ValueError, match="four columns"):
        read_probe_csv(csv)


def test_read_probe_csv_rejects_missing_magic(tmp_path):
    csv = tmp_path / "not_probe.csv"
    csv.write_text("pressure,1,0.0,1.0\n")
    with pytest.raises(ValueError, match="not a CFDX probe CSV"):
        read_probe_csv(csv)


def test_read_probe_csv_rejects_unsupported_version(tmp_path):
    csv = tmp_path / "v9.csv"
    csv.write_text(f"# {PROBE_CSV_MAGIC} v9\n# columns: probe,iteration,time,value\nx,1,0,0.0\n")
    with pytest.raises(ValueError, match="unsupported probe CSV version"):
        read_probe_csv(csv)


def test_read_probe_csv_rejects_non_numeric_value(tmp_path):
    csv = tmp_path / "nan_value.csv"
    csv.write_text(f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n# columns: probe,iteration,time,value\npressure,1,0.0,notanumber\n")
    with pytest.raises(ValueError, match="non-numeric"):
        read_probe_csv(csv)


def test_read_probe_csv_rejects_no_samples(tmp_path):
    csv = tmp_path / "empty.csv"
    csv.write_text(f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n# columns: probe,iteration,time,value\n")
    with pytest.raises(ValueError, match="no samples"):
        read_probe_csv(csv)


def test_read_probe_csv_rejects_non_monotonic_samples(tmp_path):
    csv = tmp_path / "non_mono.csv"
    csv.write_text(f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n# columns: probe,iteration,time,value\npressure,2,0.2,1.0\npressure,1,0.1,1.0\n")
    with pytest.raises(ValueError, match="monotonic"):
        read_probe_csv(csv)


@pytest.mark.parametrize("bad_name", ["", "has,comma", "has\nnewline"])
def test_write_probe_csv_rejects_bad_name(tmp_path, bad_name):
    with pytest.raises(ValueError, match="probe name must not be empty" if not bad_name else "comma or newline"):
        write_probe_csv(tmp_path / "bad.csv", [
            ProbeSeries(bad_name, (), (ProbeSample(1, 0.0, 1.0),)),
        ])


def test_write_probe_csv_rejects_empty_series(tmp_path):
    with pytest.raises(ValueError, match="at least one probe series"):
        write_probe_csv(tmp_path / "empty.csv", [])


def test_write_probe_csv_rejects_nonfinite_value(tmp_path):
    with pytest.raises(ValueError, match="non-finite"):
        write_probe_csv(tmp_path / "bad.csv", [
            ProbeSeries("p", (), (ProbeSample(1, 0.0, math.inf),)),
        ])

