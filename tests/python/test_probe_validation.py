import math
import textwrap

import pytest

from cfdx.probe_validation import (
    PROBE_CSV_MAGIC,
    PROBE_CSV_VERSION,
    ProbeSample,
    ProbeSeries,
    _format_probe_csv_value,
    _format_probe_metadata_line,
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


# Canonical probe-CSV float encoding: shortest round-trip representation with
# Python repr notation (decimal for 10^-4 <= |v| < 10^16, scientific otherwise,
# trailing ".0" for integer magnitudes). The C++ writer is implemented to emit
# these exact tokens, so this table is the shared cross-language contract between
# cfdx.io.probe_csv and cfdx.probe_validation.
CANONICAL_FLOATS = [
    (0.0, "0.0"),
    (1.0, "1.0"),
    (-0.0, "-0.0"),
    (0.1, "0.1"),
    (1.0 / 3.0, "0.3333333333333333"),
    (1e-9, "1e-09"),
    (1e-4, "0.0001"),
    (1e-5, "1e-05"),
    (1e6, "1000000.0"),
    (1e15, "1000000000000000.0"),
    (1e16, "1e+16"),
    (1e17, "1e+17"),
    (1e-308, "1e-308"),
    (123.45678901234568, "123.45678901234568"),
    (math.pi, "3.141592653589793"),
]


def test_canonical_float_format_matches_cpp():
    for value, expected in CANONICAL_FLOATS:
        assert _format_probe_csv_value(value) == expected, value


def test_write_probe_csv_value_encoding_is_canonical(tmp_path):
    series = [
        ProbeSeries("probe", (), (
            ProbeSample(i, 0.0, value) for i, (value, _) in enumerate(CANONICAL_FLOATS, start=1)
        )),
    ]
    path = write_probe_csv(tmp_path / "probes.csv", series)
    rows = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    for (value, expected), row in zip(CANONICAL_FLOATS, rows):
        # The value column must equal the canonical token shared with the C++
        # writer, and parse back to the identical double.
        assert row.split(",")[3] == expected, (value, row)
        assert float(row.split(",")[3]) == value


def test_write_probe_csv_byte_stable(tmp_path):
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
        + "".join(
            f"u_mag,{i},{_format_probe_csv_value(t)},{_format_probe_csv_value(v)}\n"
            for i, (t, v)
            in zip((1, 2, 3), ((1.0, 0.0), (1.0 / 3.0, 1.0 / 3.0), (1e-9, 123.45678901234568)))
        )
    )
    assert path.read_text() == expected


def test_write_probe_csv_overwrites_existing_file(tmp_path):
    path = tmp_path / "overwrite.csv"
    write_probe_csv(path, [
        ProbeSeries("a", (), (ProbeSample(1, 0.0, 1.0),)),
    ])
    assert path.read_text().splitlines()[-1] == "a,1,0.0,1.0"
    # Overwrite with different content; no .tmp may be left behind.
    write_probe_csv(path, [
        ProbeSeries("b", (), (ProbeSample(1, 0.0, 2.0),)),
    ])
    assert path.read_text().splitlines()[-1] == "b,1,0.0,2.0"
    assert not (tmp_path / "overwrite.csv.tmp").exists()


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


def test_read_probe_csv_rejects_bad_columns_header(tmp_path):
    csv = tmp_path / "bad_columns.csv"
    csv.write_text(f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\n# columns: foo,bar,baz,qux\nx,1,0,0.0\n")
    with pytest.raises(ValueError, match="column header mismatch"):
        read_probe_csv(csv)


def test_read_probe_csv_rejects_missing_columns_header(tmp_path):
    csv = tmp_path / "no_columns.csv"
    csv.write_text(f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}\nx,1,0.0,0.0\n")
    with pytest.raises(ValueError, match="missing the # columns: header"):
        read_probe_csv(csv)


def test_probe_metadata_line_matches_native_writer_format() -> None:
    # The C++ writer emits `# probe <name> at <x>,<y>,<z> field=<cli> unit=<u>`
    # with canonical float coordinates; the Python helper reproduces it verbatim.
    assert _format_probe_metadata_line("pressure", (0.5, 0.5, 0.5), "p", "Pa") == (
        "# probe pressure at 0.5,0.5,0.5 field=p unit=Pa"
    )
    tiny = 1e-9
    assert _format_probe_metadata_line("u", (tiny, 1e16, 0.0), "u_mag", "m/s") == (
        "# probe u at 1e-09,1e+16,0.0 field=u_mag unit=m/s"
    )


def test_write_probe_csv_emits_probe_metadata(tmp_path) -> None:
    path = write_probe_csv(tmp_path / "probes.csv", [
        ProbeSeries("pressure", (0.5, 0.5, 0.5), (
            ProbeSample(1, 0.1, 1.0),
        ), field="p", unit="Pa"),
        ProbeSeries("u_mag", (1.0, 0.0, 0.0), (
            ProbeSample(1, 0.1, 0.5),
        ), field="u_mag", unit="m/s"),
    ])
    lines = path.read_text().splitlines()
    # Metadata lines sit between the columns header and the first data row,
    # ordered by probe name to match the native writer.
    assert lines[0] == f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}"
    assert lines[1] == "# columns: probe,iteration,time,value"
    assert lines[2] == "# probe pressure at 0.5,0.5,0.5 field=p unit=Pa"
    assert lines[3] == "# probe u_mag at 1.0,0.0,0.0 field=u_mag unit=m/s"
    assert lines[4] == "pressure,1,0.1,1.0"
    assert lines[5] == "u_mag,1,0.1,0.5"


def test_read_probe_csv_round_trips_metadata(tmp_path) -> None:
    original = [
        ProbeSeries("pressure", (0.5, 0.5, 0.5), (
            ProbeSample(1, 0.1, 1.0),
            ProbeSample(2, 0.2, 1.1),
        ), field="p", unit="Pa"),
        ProbeSeries("u_mag", (1.0, 0.0, 0.0), (
            ProbeSample(1, 0.1, 0.5),
        ), field="u_mag", unit="m/s"),
    ]
    path = write_probe_csv(tmp_path / "probes.csv", original)
    loaded = read_probe_csv(path)
    assert [(s.name, s.point, s.field, s.unit) for s in loaded] == [
        ("pressure", (0.5, 0.5, 0.5), "p", "Pa"),
        ("u_mag", (1.0, 0.0, 0.0), "u_mag", "m/s"),
    ]
    for original_series, loaded_series in zip(original, loaded):
        assert loaded_series.samples == original_series.samples


def test_read_probe_csv_ignores_metadata_when_unavailable(tmp_path) -> None:
    csv = tmp_path / "legacy.csv"
    csv.write_text(textwrap.dedent("""\
        # cfdx-probe-csv v1
        # columns: probe,iteration,time,value
        pressure,1,0.5,0.3333333333333333
    """))
    loaded = read_probe_csv(csv)
    assert loaded[0].point == ()
    assert loaded[0].field is None
    assert loaded[0].unit is None

