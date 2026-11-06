"""Reproducible FGMRES/CG AMG qualification on *singular* FV Neumann matrices.

This executes the CFDX C++ test binary and parses its true-residual, error and
face-balance diagnostics. It never converts solver failure into qualification.
Times are setup + solve wall seconds, interpreted as local-machine medians only.
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import statistics
import subprocess
import sys

FIELDS = ["nx", "ny", "ratio", "direction", "krylov", "preconditioner",
          "status", "iterations", "relative_residual", "solution_error",
          "face_balance_error", "gauge_component", "seconds"]


def measure(binary: pathlib.Path, n: int, timeout: int) -> list[dict[str, str]]:
    command = [str(binary), "--benchmark", str(n), str(n)]
    run = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, timeout=timeout, check=False)
    if run.returncode != 0:
        raise RuntimeError(f"qualification binary exited {run.returncode}:\n{run.stdout[-5000:]}")
    parsed = []
    for line in run.stdout.splitlines():
        if not line.startswith("neumann_qualification_csv,"):
            continue
        values = line.split(",")[1:]
        if len(values) != len(FIELDS):
            raise RuntimeError(f"unexpected CSV schema: {line}")
        record = dict(zip(FIELDS, values))
        if int(record["nx"]) != n or int(record["ny"]) != n:
            continue  # auxiliary small-grid GMRES oracle in the same executable
        if record["status"] != "0" and record["preconditioner"] != "anisotropic_line":
            raise RuntimeError(f"AMG/CG solve failed strict gate: {line}")
        if record["status"] == "0" and float(record["relative_residual"]) > 1e-9:
            raise RuntimeError(f"uncertified convergence: {line}")
        parsed.append(record)
    if len(parsed) != 16:
        raise RuntimeError(f"expected 16 anisotropic combinations, received {len(parsed)}")
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=pathlib.Path, required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=[64, 128])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--raw", type=pathlib.Path, required=True)
    parser.add_argument("--medians", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if args.repetitions < 1 or any(n < 4 or n > 256 for n in args.sizes):
        parser.error("invalid repetitions or grid sizes")
    rows = []
    for n in args.sizes:
        for rep in range(args.repetitions):
            for record in measure(args.binary.resolve(), n, args.timeout):
                rows.append({"repetition": rep + 1, **record})
            print(f"qualified {n}x{n} repetition {rep + 1}/{args.repetitions}", flush=True)
    args.raw.parent.mkdir(parents=True, exist_ok=True)
    with args.raw.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["repetition", *FIELDS]); writer.writeheader(); writer.writerows(rows)

    groups: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in rows:
        key = tuple(row[k] for k in ("nx", "ny", "ratio", "direction", "krylov", "preconditioner"))
        groups.setdefault(key, []).append(row)
    header = ["nx", "ny", "ratio", "direction", "krylov", "preconditioner",
              "successful_runs", "total_runs", "median_iterations", "median_seconds",
              "worst_relative_residual", "worst_solution_error", "worst_face_balance_error"]
    with args.medians.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=header); writer.writeheader()
        for key, records in sorted(groups.items()):
            writer.writerow(dict(zip(header[:6], key)) | {
                "successful_runs": sum(r["status"] == "0" for r in records),
                "total_runs": len(records),
                "median_iterations": statistics.median(int(r["iterations"]) for r in records),
                "median_seconds": statistics.median(float(r["seconds"]) for r in records),
                "worst_relative_residual": max(float(r["relative_residual"]) for r in records),
                "worst_solution_error": max(float(r["solution_error"]) for r in records),
                "worst_face_balance_error": max(float(r["face_balance_error"]) for r in records),
            })
    print(f"Wrote {len(rows)} raw measurements and {len(groups)} configurations")
    return 0


if __name__ == "__main__":
    sys.exit(main())
