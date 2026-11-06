import pytest
from cfdx.probe_stream import ProbeStreamSample, parse_probe_stream_line


def test_parse_probe_stream_sample_decodes_utf8_name_and_numeric_fields() -> None:
    parsed = parse_probe_stream_line(
        "CFDX_PROBE_SAMPLE_V1\t70726f62652dceb1\tp\tPa\t7\t0.0\t101325.25"
    )

    assert parsed == ProbeStreamSample("probe-α", "p", "Pa", 7, 0.0, 101325.25)


def test_parse_probe_stream_line_ignores_regular_solver_output() -> None:
    assert parse_probe_stream_line("Iteration 7 Time = 0.0") is None
    assert parse_probe_stream_line("CFDX_PROBE_SAMPLE_V10\tanything") is None


@pytest.mark.parametrize(
    "line, message",
    [
        ("CFDX_PROBE_SAMPLE_V1", "7 columns"),
        ("CFDX_PROBE_SAMPLE_V1\tzz\tp\tPa\t1\t0\t2", "name encoding"),
        ("CFDX_PROBE_SAMPLE_V1\tff\tp\tPa\t1\t0\t2", "UTF-8"),
        ("CFDX_PROBE_SAMPLE_V1\t70\t\tPa\t1\t0\t2", "field and unit"),
        ("CFDX_PROBE_SAMPLE_V1\t70\tp\tPa\t-1\t0\t2", "negative"),
        ("CFDX_PROBE_SAMPLE_V1\t70\tp\tPa\t1\tnan\t2", "finite"),
    ],
)
def test_parse_probe_stream_line_rejects_malformed_records(
    line: str, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        parse_probe_stream_line(line)
