"""Generate a deterministic volume-skewness diagnostic JSON artifact."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "cfdx" / "python"))

import h5py
import numpy as np
from cfdx.volume_mesh_quality import analyze_cfdx_volume_skewness


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    points = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 2.0],
        ]
    )
    faces = [
        [1, 2, 3],
        [0, 2, 1],
        [0, 1, 3],
        [0, 3, 2],
        [4, 2, 1],
        [4, 1, 3],
        [4, 3, 2],
    ]
    arrays = {
        "points": points,
        "face_vertices": np.asarray(
            [v for face in faces for v in face], dtype=np.int64
        ),
        "face_offsets": np.arange(0, 24, 3, dtype=np.int64),
        "owner": np.asarray([0, 0, 0, 0, 1, 1, 1], dtype=np.int64),
        "neighbour": np.asarray([1, -1, -1, -1, -1, -1, -1], dtype=np.int64),
        "cell_faces": np.asarray([0, 1, 2, 3, 0, 4, 5, 6], dtype=np.int64),
        "cell_offsets": np.asarray([0, 4, 8], dtype=np.int64),
    }
    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / "volume-quality.h5"
        with h5py.File(source, "w") as h5:
            for name, values in arrays.items():
                h5.create_dataset(name, data=values)
        start = perf_counter()
        report = analyze_cfdx_volume_skewness(source, threshold=0.01)
        elapsed = perf_counter() - start
    if report.internal_face_count != 1 or len(report.bad_cell_ids) != 2:
        raise RuntimeError("synthetic skewness fixture did not identify both bad cells")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "fixture": "two-tetrahedra-shared-skewed-face",
                "cells": int(report.cell_skewness.size),
                "internal_faces": report.internal_face_count,
                "threshold": report.threshold,
                "max_skewness": report.max_skewness,
                "p95_skewness": report.p95_skewness,
                "bad_cell_ids": list(report.bad_cell_ids),
                "cell_neighbours": [list(values) for values in report.cell_neighbours],
                "elapsed_seconds": elapsed,
                "metric_definition": (
                    "CFDX core convention: distance from face centre to "
                    "owner-neighbour centre line / centre distance"
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
