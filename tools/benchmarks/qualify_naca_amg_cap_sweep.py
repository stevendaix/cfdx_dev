"""Fail-closed Boomer AMG interpolation-cap tuner on NACA pressure captures.

This is an offline decision aid, not an automatic CFD solver setting.  Every
candidate is solved from the same recorded initial guess and RHS and checked
against the original pressure matrix at a strict true residual tolerance.
Recommendations require a measurable median speedup; there is no tuning using
a surrogate residual or a synthetic diffusion matrix.
"""

import argparse
import csv
import io
import math
from pathlib import Path
import statistics
import subprocess
import sys


def read_replay(text: str, manifests: list[Path], tolerance: float) -> dict:
    """Reject malformed/incomplete/falsely passing replay results."""
    rows = list(csv.DictReader(io.StringIO(text)))
    if len(rows) != len(manifests):
        raise ValueError(f"expected {len(manifests)} pressure corrections, got {len(rows)}")
    metrics = {"total_ms": 0.0, "full_setup_ms": 0.0,
               "numeric_update_ms": 0.0, "krylov_ms": 0.0,
               "iterations": 0, "max_true_residual": 0.0}
    for index, (record, manifest) in enumerate(zip(rows, manifests), start=1):
        if Path(record.get("manifest", "")).resolve() != manifest.resolve():
            raise ValueError(f"pressure capture order/identity mismatch at {index}")
        if record.get("status") != "PASS" or record.get("sequence") != str(index):
            raise ValueError(f"not a certified pressure solve at {index}")
        true_residual = float(record["true_relative_residual"])
        if not math.isfinite(true_residual) or true_residual < 0 or true_residual > tolerance:
            raise ValueError(f"true residual not satisfied at {index}: {true_residual}")
        if int(record["rows"]) <= 0 or int(record["iterations"]) < 0:
            raise ValueError("invalid linear system/iteration count")
        for field in ("total_ms", "full_setup_ms", "numeric_update_ms", "krylov_ms"):
            value = float(record[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid {field} measurement")
            metrics[field] += value
        metrics["iterations"] += int(record["iterations"])
        metrics["max_true_residual"] = max(metrics["max_true_residual"], true_residual)
    metrics["operator_complexity"] = float(rows[-1]["amg_operator_complexity"])
    metrics["levels"] = int(rows[-1]["amg_levels"])
    if not math.isfinite(metrics["operator_complexity"]) or metrics["operator_complexity"] < 1:
        raise ValueError("invalid operator complexity")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--captures", type=Path, required=True)
    parser.add_argument("--caps", default="0,3,4,5,6", help="0=unlimited; each other value in [2,16]")
    parser.add_argument("--repetitions", type=int, default=4)
    parser.add_argument("--min-speedup", type=float, default=1.10)
    parser.add_argument("--true-tolerance", type=float, default=1e-10)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    caps = [int(item) for item in args.caps.split(",")]
    if len(set(caps)) != len(caps) or 0 not in caps or any(cap != 0 and not 2 <= cap <= 16 for cap in caps):
        parser.error("caps must be distinct, include 0, and otherwise be in [2,16]")
    if not 3 <= args.repetitions <= 20 or not 1.0 < args.min_speedup < 10.0:
        parser.error("require 3..20 repetitions and min-speedup in (1,10)")
    if not 0 < args.true_tolerance <= 1e-9 or args.timeout <= 0:
        parser.error("unsafe residual threshold or timeout")
    if not args.replay.is_file():
        parser.error("replay binary not found")

    raw = []
    summaries = []
    for mesh in ("N128", "N256"):
        captures = sorted(args.captures.glob(f"{mesh}_portable*.manifest"))
        if len(captures) != 3:
            raise ValueError(f"expected three captured pressure corrections for {mesh}")
        for repetition in range(args.repetitions):
            # Rotate cap order: baseline must not systematically profit from
            # cold caches or from always running first or last.
            order = caps[repetition % len(caps):] + caps[:repetition % len(caps)]
            for cap in order:
                method = "fgmres_boomer_fast" + (f"_cap{cap}" if cap else "")
                command = [str(args.replay.resolve()), method, "reuse", "1",
                           *(str(c.resolve()) for c in captures)]
                result = subprocess.run(command, capture_output=True,
                                        text=True, check=False, timeout=args.timeout)
                if result.returncode != 0:
                    raise RuntimeError(f"{mesh}/{method}: replay exited {result.returncode}: {result.stderr[-800:]}")
                measurements = read_replay(result.stdout, captures, args.true_tolerance)
                raw.append(dict(mesh=mesh, cap=cap, repetition=repetition, **measurements))
                print(f"{mesh} cap={cap} rep={repetition} "
                      f"time={measurements['total_ms']:.2f} ms "
                      f"max_true_residual={measurements['max_true_residual']:.3g}", flush=True)
        median_by_cap = {}
        for cap in caps:
            records = [r for r in raw if r["mesh"] == mesh and r["cap"] == cap]
            entry = dict(mesh=mesh, cap=cap, samples=len(records))
            for key in ("total_ms", "full_setup_ms", "numeric_update_ms", "krylov_ms", "operator_complexity"):
                entry[key] = statistics.median(float(r[key]) for r in records)
            entry["iterations"] = int(statistics.median(r["iterations"] for r in records))
            entry["max_true_residual"] = max(r["max_true_residual"] for r in records)
            median_by_cap[cap] = entry
        baseline = median_by_cap[0]["total_ms"]
        fastest = min(caps, key=lambda cap: median_by_cap[cap]["total_ms"])
        # Do not recommend a new cap just because of noise: require a minimum
        # observed end-to-end pressure-replay improvement over the baseline.
        preferred = fastest if baseline >= args.min_speedup * median_by_cap[fastest]["total_ms"] else 0
        for cap in caps:
            entry = median_by_cap[cap]
            entry["speedup_vs_default"] = baseline / entry["total_ms"]
            entry["offline_recommended"] = cap == preferred
            entry["qualification"] = "PASS_LINEAR_TRUE_RESIDUAL_ONLY"
            summaries.append(entry)
        print(f"{mesh}: offline recommended cap={preferred} "
              f"(unlimited=0), candidate cap={fastest}, "
              f"gain={baseline/median_by_cap[fastest]['total_ms']:.3f}x", flush=True)
    for path, records in ((args.raw, raw), (args.summary, summaries)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    print("PASS: all candidates certified; offline recommendations only", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, KeyError, OSError,
            subprocess.TimeoutExpired) as exc:
        print(f"AMG_CAP_QUALIFICATION_FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
