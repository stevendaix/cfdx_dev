"""Benchmark metadata-only inspection and bounded samples on sparse HDF5 meshes.

The generated datasets have production-scale logical dimensions but contain no
materialized mesh payload. This isolates HDF5 metadata traversal and bounded
sampling; it is not a PyVista rendering benchmark.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from time import perf_counter
from typing import Any

import h5py

from cfdx.case_inspector import inspect_cfdx_file, read_hdf5_dataset_sample


def _peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    except (ImportError, AttributeError):
        return None
    return float(peak / (1024 * 1024 if sys.platform == "darwin" else 1024))


def _create_sparse_case(path: Path, n_cells: int) -> None:
    n_points = 4 * n_cells
    n_faces = 6 * n_cells
    n_connectivity = 4 * n_faces
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["n_cells"] = n_cells
        h5.attrs["source_case_name"] = f"sparse-{n_cells}-cells"
        h5.create_dataset(
            "points",
            shape=(n_points, 3),
            dtype="f8",
            chunks=(1024, 3),
            fill_time="never",
        )
        for name, dtype in (
            ("owner", "i8"),
            ("neighbour", "i8"),
            ("patch_face_ids", "u8"),
        ):
            h5.create_dataset(
                name,
                shape=(n_faces,),
                dtype=dtype,
                chunks=(4096,),
                fill_time="never",
            )
        h5.create_dataset(
            "face_vertices",
            shape=(n_connectivity,),
            dtype="u8",
            chunks=(4096,),
            fill_time="never",
        )
        h5.create_dataset(
            "face_offsets",
            shape=(n_faces + 1,),
            dtype="u8",
            chunks=(4096,),
            fill_time="never",
        )
        h5.create_dataset("patch_face_offsets", data=[0, n_faces])
        h5.create_dataset(
            "fields/scalar/p",
            shape=(n_cells,),
            dtype="f8",
            chunks=(4096,),
            fill_time="never",
        )


def _benchmark_one(n_cells: int) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="cfdx-hdf5-bench-") as directory:
        path = Path(directory) / f"case-{n_cells}.cfdx.h5"
        _create_sparse_case(path, n_cells)
        rss_before = _peak_rss_mb()

        start = perf_counter()
        summary = inspect_cfdx_file(path)
        inspect_ms = (perf_counter() - start) * 1000

        start = perf_counter()
        sample = read_hdf5_dataset_sample(path, "/points", max_values=32)
        sample_ms = (perf_counter() - start) * 1000
        rss_after = _peak_rss_mb()

        if summary.n_cells != n_cells or summary.preview_supported:
            raise RuntimeError("large sparse fixture did not hit metadata-only mode")
        if not sample.values or len(sample.values) > 32:
            raise RuntimeError("bounded dataset sample violated its 32-value cap")

        logical_bytes = sum(item.size_bytes or 0 for item in summary.datasets)
        return {
            "cells": n_cells,
            "points": summary.n_points,
            "faces": summary.n_faces,
            "logical_dataset_bytes": logical_bytes,
            "physical_file_bytes": path.stat().st_size,
            "inspect_ms": round(inspect_ms, 3),
            "sample_ms": round(sample_ms, 3),
            "sample_values": len(sample.values),
            "peak_rss_mb": round(rss_after, 3) if rss_after is not None else None,
            "rss_growth_mb": (
                round(max(0.0, rss_after - rss_before), 3)
                if rss_before is not None and rss_after is not None
                else None
            ),
            "preview_block_reason": summary.preview_block_reason,
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sizes",
        default="1000000,10000000,50000000",
        help="comma-separated logical cell counts",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker is not None:
        print(json.dumps(_benchmark_one(args.worker)))
        return 0

    sizes = [int(value) for value in args.sizes.split(",") if value.strip()]
    if not sizes or any(size < 1 for size in sizes):
        parser.error("--sizes must contain positive integer cell counts")

    results = []
    for size in sizes:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker", str(size)],
            check=True,
            capture_output=True,
            text=True,
            env=os.environ.copy(),
        )
        results.append(json.loads(completed.stdout))

    report = {
        "scope": "metadata inspection and bounded HDF5 sampling only; no rendering",
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "machine": platform.machine(),
            "logical_cpus": os.cpu_count(),
            "h5py": h5py.__version__,
            "hdf5": h5py.version.hdf5_version,
        },
        "results": results,
    }
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
