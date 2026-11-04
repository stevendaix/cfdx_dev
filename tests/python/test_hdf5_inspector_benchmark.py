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
    assert result["sample_values"] <= 32
    assert result["preview_block_reason"]
    assert result["physical_file_bytes"] < result["logical_dataset_bytes"]
