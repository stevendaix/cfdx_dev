"""Batch-compare CFDX preconditioners on EXACT captured original CSR systems.

No candidate is silently discarded: missing files, nonconverged solves, bad
CSV or nonfinite timings make the entire qualification return nonzero.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import subprocess
import sys
from pathlib import Path

DEFAULT_METHODS = ("jacobi", "ilu0", "native", "directional", "energy")


def compare(
    directory: Path, binary: Path, methods: tuple[str, ...], timeout: float
) -> tuple[list[dict[str, str]], list[str]]:
    files = sorted(directory.glob("*.cfdxlin"))
    if not files:
        return [], ["no captured systems found"]
    rows: list[dict[str, str]] = []
    errors: list[str] = []
    for path in files:
        try:
            result = subprocess.run(
                [str(binary.resolve()), str(path.resolve()), *methods],
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
            parsed = list(csv.DictReader(io.StringIO(result.stdout)))
            reasons: dict[str, str] = {}
            for line in result.stderr.splitlines():
                for prefix in ("PRECONDITIONER_SETUP_ERROR ", "REPLAY_SOLVE_FAILED "):
                    if line.startswith(prefix):
                        method, separator, explanation = line[len(prefix) :].partition(
                            ": "
                        )
                        if separator and method in methods:
                            reasons[method] = explanation.strip()
            names = [row.get("preconditioner", "") for row in parsed]
            if names != list(methods):
                raise ValueError("method rows do not match the requested list")
            if result.returncode != 0:
                errors.append(f"{path.name}: replay returned {result.returncode}")
            for item in parsed:
                if item.get("status") != "converged":
                    errors.append(
                        f"{path.name}: {item.get('preconditioner')} not converged"
                    )
                try:
                    elapsed = float(item["total_ms"])
                    residual = float(item["true_relative"])
                    iterations = int(item["iterations"])
                except (TypeError, ValueError, KeyError) as exc:
                    raise ValueError("bad numeric replay evidence") from exc
                if not all(map(math.isfinite, (elapsed, residual))) or (
                    elapsed < 0 or residual < 0 or iterations < 0
                ):
                    errors.append(f"{path.name}: nonfinite or negative replay evidence")
                method = item.get("preconditioner", "")
                reason = reasons.get(method, "")
                if item.get("status") != "converged" and not reason:
                    reason = f"status={item.get('status', 'unknown')}"
                rows.append({"capture": path.name, **item, "failure_reason": reason})
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            errors.append(f"{path.name}: {exc}")
    return rows, errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--method", action="append", dest="methods")
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()
    methods = tuple(args.methods or DEFAULT_METHODS)
    if (
        len(methods) != len(set(methods))
        or not methods
        or (not math.isfinite(args.timeout) or args.timeout <= 0)
    ):
        parser.error("invalid method selection or timeout")
    rows, errors = compare(args.directory, args.binary, methods, args.timeout)
    header = [
        "capture",
        "preconditioner",
        "krylov",
        "status",
        "iterations",
        "total_ms",
        "true_abs_l2",
        "true_relative",
        "failure_reason",
    ]
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    for error in errors:
        print("CFDX_REPLAY_BATCH_FAIL:", error, file=sys.stderr)
    print(
        f"captures={len({r['capture'] for r in rows})} "
        f"candidate_results={len(rows)} failures={len(errors)}"
    )
    return 0 if rows and not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
