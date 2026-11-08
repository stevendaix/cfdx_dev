"""Compare native PCD Fp preconditioners on exactly the SAME captured systems.

Each capture is a byte-exact production CSR, RHS and initial iterate. Results
are diagnostic and never qualify the COUPLED solver or a physical benchmark.
No silent fallbacks, fictitious convergence, or changed physical thresholds.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

METHODS = (
    "ilu0",
    "ilu0_equilibrated",
    "amg",
    "amg_robust",
    "amg_equilibrated",
    "amg_sa",
    "amg_lines",
    "amg_semicoarsen",
    "jacobi",
)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as src:
        while block := src.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def compare(
    captures: list[Path],
    executable: Path,
    *,
    budget: int = 100,
    tolerance: float = 1e-8,
    runner=subprocess.run,
) -> dict:
    if not executable.is_file() or not captures:
        raise ValueError("missing replay executable or production captures")
    if budget < 1 or not (1e-14 <= tolerance <= 1e-4):
        raise ValueError("invalid replay controls")
    records = []
    for capture in captures:
        if not capture.is_file():
            raise FileNotFoundError(capture)
        invocation = runner(
            [str(executable), str(capture), str(budget), str(tolerance)],
            text=True,
            capture_output=True,
            check=False,
        )
        if invocation.returncode:
            raise RuntimeError(f"{capture}: replay failed: {invocation.stderr}")
        rows = [
            json.loads(line) for line in invocation.stdout.splitlines() if line.strip()
        ]
        if len(rows) != len(METHODS) or {r["method"] for r in rows} != set(METHODS):
            raise ValueError(f"{capture}: incomplete or duplicate method report")
        for row in rows:
            if row["rows"] <= 0 or row["nnz"] < row["rows"]:
                raise ValueError(f"{capture}: invalid Fp dimensions")
            if row["method"] not in METHODS:
                raise ValueError(f"{capture}: unknown method")
            record = {**row, "capture_sha256": digest(capture)}
            record["first_meets_true_tolerance"] = (
                row["setup_ok"]
                and row["first_status"] == 0
                and row["first_true_relative_rhs"] <= tolerance
            )
            record["two_pass_meets_true_tolerance"] = (
                row["setup_ok"]
                and row["best_true_relative_rhs"] <= tolerance
                and (
                    row["first_status"] == 0
                    or (
                        row["second_attempted"]
                        and row["second_improved"]
                        and row["second_status"] == 0
                    )
                )
            )
            records.append(record)
    return {
        "schema": "cfdx-pcd-fp-replay/v1",
        "method": "original_CSR_and_rhs_independent_true_residual",
        "replay_executable_sha256": digest(executable),
        "inner_max_iterations_per_pass": budget,
        "inner_tolerance": tolerance,
        "capture_count": len(captures),
        "records": records,
        "note": (
            "Results qualify the Fp replay only; no Schur or physical NACA/N9"
            " qualification, no convergence inferred from CTest or history."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--capture", type=Path, action="append", default=[])
    ap.add_argument("--capture-dir", type=Path)
    ap.add_argument("--executable", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--budget", type=int, default=100)
    ap.add_argument("--tolerance", type=float, default=1e-8)
    ap.add_argument("--require-all-true-residuals", action="store_true")
    args = ap.parse_args(argv)
    try:
        captures = list(args.capture)
        if args.capture_dir:
            captures += sorted(args.capture_dir.glob("*.cfp"))
        if len(captures) != len(set(captures)):
            raise ValueError("duplicate capture paths")
        report = compare(
            captures, args.executable, budget=args.budget, tolerance=args.tolerance
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "pcd_fp_replay.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n"
        )
        cols = [
            "capture",
            "capture_sha256",
            "rows",
            "nnz",
            "method",
            "setup_ok",
            "first_iterations",
            "first_true_relative_rhs",
            "best_true_relative_rhs",
            "second_attempted",
            "second_improved",
            "first_meets_true_tolerance",
            "two_pass_meets_true_tolerance",
            "total_seconds_two_passes",
            "probe_A_correction_l2",
            "probe_cosine",
            "bounded_action_fallbacks",
            "error",
        ]
        with (args.output_dir / "pcd_fp_replay.csv").open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(report["records"])
        passes = sum(r["two_pass_meets_true_tolerance"] for r in report["records"])
        print(
            f"PCD_FP_REPLAY: captures={len(captures)} entries={len(report['records'])} "
            f"strict_inner_pass={passes} (NOT physical certification)"
        )
        if args.require_all_true_residuals and passes != len(report["records"]):
            return 2
        return 0
    except (ValueError, OSError, RuntimeError, KeyError, json.JSONDecodeError) as exc:
        print(f"PCD_FP_REPLAY: ERROR {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
