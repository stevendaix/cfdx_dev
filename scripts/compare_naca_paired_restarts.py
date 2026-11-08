"""Fail-closed comparison of two same-mesh, accepted NACA U,p,phi snapshots."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def fnv64(path: Path) -> str:
    value = 14695981039346656037
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            for byte in block:
                value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return f"{value:016x}"


def load_manifest(path: Path) -> dict:
    lines = path.read_text().splitlines()
    if len(lines) != 2 or lines[0] != "CFDX_NACA_PAIRED_V2":
        raise ValueError("invalid NACA paired-checkpoint manifest")
    fields = lines[1].split(",")
    if len(fields) != 9:
        raise ValueError("invalid paired manifest column count")
    step, cells, faces = (int(t) for t in fields[:3])
    momentum_reference, monitor_reference = map(float, fields[3:5])
    if step <= 0 or cells <= 0 or faces <= 0 or not all(
        math.isfinite(x) and x > 0 for x in (momentum_reference, monitor_reference)
    ):
        raise ValueError("invalid iteration or monitor references")
    for filename in fields[5:7]:
        if not filename or ".." in filename or Path(filename).name != filename:
            raise ValueError("unsafe payload path")
    paths = [path.parent / fields[5], path.parent / fields[6]]
    for fp, expected in zip(paths, fields[7:]):
        if fnv64(fp) != expected:
            raise ValueError(f"checkpoint checksum mismatch: {fp}")
    return dict(iteration=step, n_cells=cells, n_faces=faces,
                momentum_reference=momentum_reference,
                nonlinear_reference=monitor_reference, files=paths)


def compare(a: Path, b: Path, *, tolerance: float = 1e-12) -> dict:
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("invalid strict comparison tolerance")
    left, right = load_manifest(a), load_manifest(b)
    if any(left[k] != right[k] for k in (
            "iteration", "n_cells", "n_faces", "momentum_reference",
            "nonlinear_reference")):
        raise ValueError(
            "different iteration, mesh count, or convergence normalization")
    results = {}
    for kind, fieldnames, coordinate_columns, n_expected, ix in (
        ("cell", ("ux", "uy", "p"), ("cell", "x", "y", "z"),
         left["n_cells"], 0),
        ("face", ("mass_flux",),
         ("face", "owner", "neighbour", "cx", "cy", "cz", "sx", "sy", "sz"),
         left["n_faces"], 1),
    ):
        with (left["files"][ix].open(newline="") as fa,
              right["files"][ix].open(newline="") as fb):
            rows_a, rows_b = csv.DictReader(fa), csv.DictReader(fb)
            if rows_a.fieldnames != rows_b.fieldnames or any(
                column not in rows_a.fieldnames
                for column in (*coordinate_columns, *fieldnames)
            ):
                raise ValueError("different CSV schemas")
            max_error = {key: 0.0 for key in fieldnames}
            sum_sq = {key: 0.0 for key in fieldnames}
            n = 0
            for row_a, row_b in zip(rows_a, rows_b):
                for key in coordinate_columns:
                    if row_a[key] != row_b[key]:
                        raise ValueError(
                            f"{kind} connectivity or geometry mismatch, "
                            f"row {n}, {key}")
                for key in fieldnames:
                    values = float(row_a[key]), float(row_b[key])
                    if not all(math.isfinite(v) for v in values):
                        raise ValueError(
                            "non-finite field in restart comparison")
                    error = abs(values[0] - values[1])
                    max_error[key] = max(max_error[key], error)
                    sum_sq[key] += error * error
                n += 1
            if (n != n_expected or next(rows_a, None) is not None
                    or next(rows_b, None) is not None):
                raise ValueError(
                    f"{kind} row-count mismatch: read {n}, "
                    f"expected {n_expected}")
        results[kind] = {
            key: {"linf": max_error[key],
                  "l2_rms": math.sqrt(sum_sq[key] / n)}
            for key in fieldnames
        }
    max_difference = max(
        v["linf"] for family in results.values() for v in family.values())
    return dict(iteration=left["iteration"],
                max_abs_difference=max_difference, tolerance=tolerance,
                verdict="PASS" if max_difference <= tolerance else "FAIL",
                norms=results)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--continuous", required=True, type=Path)
    parser.add_argument("--resumed", required=True, type=Path)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()
    try:
        report = compare(args.continuous, args.resumed, tolerance=args.atol)
    except (OSError, ValueError, KeyError) as exc:
        report = dict(verdict="ERROR", error=str(exc))
    print(json.dumps(report, indent=2))
    return 0 if report.get("verdict") == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
