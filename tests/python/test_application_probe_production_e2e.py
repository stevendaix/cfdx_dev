"""Application-routed probe qualification against a real production solver run.

The shared Application owns the probe CSV contract (``read_probe_results`` /
``export_probe_results``) over the strict native schema. This test qualifies
that contract against an artifact actually emitted by the production solver,
not against a synthetic ``ProbeSeries``: it runs the real solver with declared
probes, reads the emitted CSV back through the Application, and asserts the
round-trip preserves the schema and the samples the solver produced.

The production solver and mesh are provided by ctest via environment variables,
so the test skips cleanly when the fixtures are unavailable.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest
from cfdx.probe import Probe, ProbeCatalog, ProbeField
from cfdx.probe_validation import ProbeSample

from cfdx.application import Application


def _production_env() -> dict[str, str]:
    return {
        "CFDX_PRODUCTION_SOLVER": os.environ.get("CFDX_PRODUCTION_SOLVER", ""),
        "CFDX_PRODUCTION_MESH": os.environ.get("CFDX_PRODUCTION_MESH", ""),
    }


def test_application_probe_csv_round_trips_real_solver_output(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("CFDX_PRODUCTION_SOLVER/CFDX_PRODUCTION_MESH are provided by ctest")
    assert Path(solver).is_file(), f"production solver missing: {solver}"
    assert Path(mesh).is_file(), f"production mesh missing: {mesh}"

    output_dir = tmp_path / "probe_run"
    output_dir.mkdir()
    csv_path = output_dir / "case.probes.csv"

    catalog = ProbeCatalog()
    catalog.add(Probe(name="pressure", x=0.05, y=0.5, z=0.5, field=ProbeField.PRESSURE))
    catalog.add(Probe(name="ux", x=0.05, y=0.5, z=0.5, field=ProbeField.U_X))

    completed = subprocess.run(
        [
            solver,
            "--mesh",
            mesh,
            "--output-dir",
            str(output_dir),
            "--iterations",
            "5",
            "--probe",
            catalog.specs()[0],
            "--probe",
            catalog.specs()[1],
            "--probe-csv",
            str(csv_path),
        ],
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr[-4000:]
    assert csv_path.is_file(), "production solver did not emit a probe CSV"

    application = Application()
    series = application.read_probe_results(csv_path)
    assert len(series) == 2
    by_name = {probe.name: probe for probe in series}
    assert by_name["pressure"].field == "p"
    assert by_name["pressure"].unit == "Pa"
    assert by_name["ux"].field == "u_x"
    assert by_name["ux"].unit == "m/s"

    # The solver emits one sample per iteration per probe; the Application
    # reader must preserve every sample the solver wrote.
    assert all(len(probe.samples) >= 1 for probe in series)
    for probe in series:
        for sample in probe.samples:
            assert isinstance(sample, ProbeSample)
            assert sample.iteration >= 0

    # The contract is symmetric: the same schema the solver wrote must be
    # reproducible verbatim through the Application's writer.
    round_trip = tmp_path / "round_trip.probes.csv"
    application.export_probe_results(round_trip, series)
    assert round_trip.read_text(encoding="utf-8") == csv_path.read_text(
        encoding="utf-8"
    )


def test_application_probe_csv_rejects_solver_emitted_ambiguous_history(
    tmp_path: Path,
) -> None:
    """A CSV with duplicate iteration IDs for one probe must fail to read.

    The production solver writes one sample per iteration per probe, so a file
    that violates that contract is not a valid probe history and the
    Application reader must reject it rather than silently merging samples.
    """
    bad_csv = tmp_path / "ambiguous.probes.csv"
    bad_csv.write_text(
        "# cfdx-probe-csv v1\n"
        "# columns: probe,iteration,time,value\n"
        "# probe p at 0.05,0.5,0.5 field=p unit=Pa\n"
        "p,1,0.0,101325.0\n"
        "p,1,0.0,101326.0\n",
        encoding="utf-8",
    )

    application = Application()
    with pytest.raises(ValueError, match="duplicate iteration"):
        application.read_probe_results(bad_csv)
