"""Generate a reproducible HDF5 surface-quality diagnostic artifact."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from time import perf_counter

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PACKAGE_ROOT = REPOSITORY_ROOT / "src" / "cfdx" / "python"
sys.path.insert(0, str(PYTHON_PACKAGE_ROOT))

import h5py
import numpy as np

from cfdx.mesh_quality import analyze_cfdx_surface_file


def run_benchmark(output: Path, resolution: int = 128) -> dict[str, object]:
    """Build a deterministic warped quad surface and export measured diagnostics."""
    if resolution < 2 or resolution > 512:
        raise ValueError("resolution must be between 2 and 512")
    coordinates = np.linspace(0.0, 1.0, resolution + 1)
    xx, yy = np.meshgrid(coordinates, coordinates, indexing="xy")
    zz = 0.01 * np.sin(np.pi * xx) * np.sin(np.pi * yy)
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

    i, j = np.meshgrid(
        np.arange(resolution, dtype=np.int64),
        np.arange(resolution, dtype=np.int64),
        indexing="xy",
    )
    lower_left = (j * (resolution + 1) + i).ravel()
    face_vertices = np.column_stack(
        (
            lower_left,
            lower_left + 1,
            lower_left + resolution + 2,
            lower_left + resolution + 1,
        )
    ).ravel()
    face_offsets = np.arange(0, face_vertices.size + 1, 4, dtype=np.int64)

    with tempfile.TemporaryDirectory(prefix="cfdx-surface-quality-") as directory:
        path = Path(directory) / "warped-surface.cfdx.h5"
        with h5py.File(path, "w") as h5:
            h5.create_dataset("points", data=points, chunks=(4096, 3))
            h5.create_dataset("face_vertices", data=face_vertices, chunks=(16384,))
            h5.create_dataset("face_offsets", data=face_offsets, chunks=(4096,))
        start = perf_counter()
        report = analyze_cfdx_surface_file(path)
        elapsed_ms = (perf_counter() - start) * 1000.0

    expected_faces = resolution * resolution
    if (
        report.face_count != expected_faces
        or report.valid_face_count != expected_faces
    ):
        raise RuntimeError(
            "surface quality analysis did not cover every generated face"
        )
    if report.warpage_ratio_max <= 0.0:
        raise RuntimeError("warped-surface probe did not produce nonzero warpage")

    result: dict[str, object] = {
        "benchmark": "cfdx_surface_quality",
        "fixture": "synthetic_warped_quad_surface",
        "resolution": resolution,
        "face_count": report.face_count,
        "valid_face_count": report.valid_face_count,
        "degenerate_face_count": len(report.degenerate_face_ids),
        "zero_length_edge_face_count": len(report.zero_length_edge_face_ids),
        "face_area_min": report.face_area_min,
        "face_area_max": report.face_area_max,
        "edge_length_ratio_min": report.edge_length_ratio_min,
        "edge_length_ratio_p95": report.edge_length_ratio_p95,
        "edge_length_ratio_max": report.edge_length_ratio_max,
        "warpage_ratio_p95": report.warpage_ratio_p95,
        "warpage_ratio_max": report.warpage_ratio_max,
        "worst_warpage_face_ids": list(report.worst_warpage_face_ids),
        "elapsed_ms": round(elapsed_ms, 3),
        "notes": (
            "Synthetic surface-only diagnostic; does not measure volume-cell "
            "quality, mesh rendering, solver throughput or mesh repair."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/cfdx-surface-quality.json"),
    )
    parser.add_argument("--resolution", type=int, default=128)
    args = parser.parse_args()
    result = run_benchmark(args.output, args.resolution)
    print(
        f"Surface quality export: faces={result['face_count']}, "
        f"warpage_max={result['warpage_ratio_max']:.6g}, "
        f"elapsed_ms={result['elapsed_ms']:.3f}, output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
