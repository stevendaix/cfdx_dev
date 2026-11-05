"""Fail-closed, paired V4/V5 AMG replay on identical production NACA captures.

Both executables see the same ordered manifest paths. The replay results must
identify all three inputs and individually satisfy the original true residual.
An AB/BA run order limits (but does not eliminate) temperature/cache bias.
This measures linear replay time, not full CFD time-to-convergence.
"""

import argparse
import csv
import io
import math
from pathlib import Path
import statistics
import subprocess
import sys

from qualify_naca_amg_cap_sweep import read_replay

METHODS = (
    "fgmres_boomer_fast",
    "fgmres_boomer_fast_cap4",
    "fgmres_sa_fast",
)


def read_certified_replay(output: str, captures: list[Path], tolerance: float) -> dict:
    """Require identity/order, finite metrics and true residuals on every solve."""
    return read_replay(output, captures, tolerance)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--captures", required=True, type=Path)
    parser.add_argument("--repetitions", type=int, default=4)
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--true-tolerance", type=float, default=1e-10)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    methods = args.methods.split(",")
    if not 2 <= args.repetitions <= 20:
        parser.error("repetitions must be in [2, 20]")
    if not methods or len(methods) != len(set(methods)) or any(m not in METHODS for m in methods):
        parser.error("methods must be distinct supported replay methods")
    if not 0 < args.true_tolerance <= 1e-9 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("unsafe tolerance or timeout")
    if not args.base.is_file() or not args.candidate.is_file():
        parser.error("V4 and V5 replay executables must exist")

    records: list[dict] = []
    for mesh in ("N128", "N256"):
        captures = sorted(args.captures.glob(f"{mesh}_portable*.manifest"))
        if len(captures) != 3:
            raise ValueError(f"{mesh}: expected exactly 3 NACA pressure captures")
        captures = [path.resolve() for path in captures]
        for method in methods:
            for repetition in range(args.repetitions):
                versions = (("V4", args.base), ("V5", args.candidate))
                if repetition % 2:
                    versions = tuple(reversed(versions))
                for version, executable in versions:
                    command = [str(executable.resolve()), method, "reuse", "1",
                               *(str(capture) for capture in captures)]
                    process = subprocess.run(
                        command, check=False, capture_output=True,
                        text=True, timeout=args.timeout,
                    )
                    if process.returncode != 0:
                        raise RuntimeError(
                            f"{mesh}/{method}/{version}: exit {process.returncode}: "
                            f"{process.stderr[-500:]}"
                        )
                    metrics = read_certified_replay(
                        process.stdout, captures, args.true_tolerance,
                    )
                    record = dict(mesh=mesh, method=method, repeat=repetition,
                                  version=version, total_ms=metrics["total_ms"],
                                  full_setup_ms=metrics["full_setup_ms"],
                                  numeric_update_ms=metrics["numeric_update_ms"],
                                  krylov_ms=metrics["krylov_ms"],
                                  iterations_total=metrics["iterations"],
                                  max_true_residual=metrics["max_true_residual"],
                                  operator_complexity=metrics["operator_complexity"],
                                  levels=metrics["levels"])
                    records.append(record)
                    print(mesh, method, repetition, version,
                          f"{record['total_ms']:.2f} ms", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print("--- paired results ---")
    for mesh in ("N128", "N256"):
        for method in methods:
            selected = [r for r in records if r["mesh"] == mesh and r["method"] == method]
            base = statistics.median(r["total_ms"] for r in selected if r["version"] == "V4")
            opt = statistics.median(r["total_ms"] for r in selected if r["version"] == "V5")
            print(f"{mesh} {method} V4={base:.2f}ms V5={opt:.2f}ms "
                  f"speedup={base / opt:.3f}x")
    print(f"PASS: all replay rows matched captures and true residual <= {args.true_tolerance}")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, KeyError, OSError,
            subprocess.TimeoutExpired) as exc:
        print("PAIRED_BENCHMARK_FAIL:", exc, file=sys.stderr)
        raise SystemExit(1) from exc
