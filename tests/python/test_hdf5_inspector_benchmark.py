"""CLI regression tests for the sparse HDF5 inspector benchmark."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


def test_benchmark_cli_runs_without_pythonpath_or_scripts_shadowing(tmp_path):
    repository = Path(__file__).resolve().parents[2]
    script = repository / "scripts" / "benchmark_cfdx_hdf5_inspector.py"
    output = tmp_path / "benchmark.json"
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--sizes",
            "1000000",
            "--output",
            str(output),
        ],
        cwd=tmp_path,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )

    report = json.loads(output.read_text(encoding="utf-8"))
    assert json.loads(completed.stdout) == report
    assert "no rendering" in report["scope"]
    assert len(report["results"]) == 1
    result = report["results"][0]
    assert result["cells"] == 1_000_000
    assert result["sample_values"] == 30
    field_statistics = result["field_statistics"]
    assert field_statistics["path"] == "/fields/scalar/temperature_probe"
    assert field_statistics["value_count"] == 4096
    assert field_statistics["finite_count"] == 4096
    assert field_statistics["non_finite_count"] == 0
    assert field_statistics["minimum"] == -2048.0
    assert field_statistics["maximum"] == 2047.0
    assert field_statistics["scanned_bytes"] == 4096 * 8
    assert field_statistics["elapsed_ms"] >= 0
    assert result["points_allocated_bytes"] > 0
    assert result["allocated_dataset_bytes"] >= result["points_allocated_bytes"]
    storage_audit = {item["path"]: item for item in result["dataset_storage_audit"]}
    assert (
        storage_audit["/points"]["allocated_bytes"] == result["points_allocated_bytes"]
    )
    assert storage_audit["/points"]["allocated_ratio"] is not None
    assert result["preview_block_reason"]
    assert result["physical_file_bytes"] < result["logical_dataset_bytes"]
