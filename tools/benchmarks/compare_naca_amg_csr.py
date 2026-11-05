"""Paired, fail-closed NACA pressure snapshots benchmark of AMG Galerkin CSR.

Reference is the V3 binary (std::map Galerkin), and candidate is this V4
binary (symbolic/numeric CSR), with optional native Boomer CF cap=4.
No NACA forces/physical convergence claims are made by this linear replay.
"""
import argparse
import csv
import io
import math
from pathlib import Path
import statistics
import subprocess
import sys


def run(binary: Path, method: str, captures: list[Path], timeout: float) -> list[dict]:
    proc = subprocess.run(
        [str(binary.resolve()), method, "reuse", "1", *map(str, captures)],
        capture_output=True, text=True, check=False, timeout=timeout,
    )
    if proc.returncode:
        raise RuntimeError(f"{binary.name}/{method}: exit={proc.returncode} {proc.stderr[-900:]}")
    records = list(csv.DictReader(io.StringIO(proc.stdout)))
    if len(records) != len(captures):
        raise RuntimeError(f"expected {len(captures)} rows, got {len(records)}: {binary}/{method}")
    for record in records:
        residual = float(record["true_relative_residual"])
        if record["status"] != "PASS" or not math.isfinite(residual) or residual > 1e-10:
            raise RuntimeError(f"uncertified true residual: {binary}/{method}: {record}")
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--captures", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--timeout-seconds", type=float, default=90)
    args = parser.parse_args()
    if args.repetitions < 2 or args.repetitions > 20:
        parser.error("repetitions outside [2,20]")

    results = []
    summary = []
    for mesh in ("N128", "N256"):
        captures = sorted(args.captures.glob(f"{mesh}_portable*.manifest"))
        if len(captures) != 3:
            raise RuntimeError(f"expected 3 NACA {mesh} captures, got {len(captures)}")
        for method in ("fgmres_boomer_fast", "fgmres_sa_fast"):
            configs = [("v3_map", args.baseline, method),
                       ("v4_csr", args.candidate, method)]
            if method == "fgmres_boomer_fast":
                configs.append(("v4_csr_cap4", args.candidate,
                                "fgmres_boomer_fast_cap4"))
            grouped = {key: [] for key, _, _ in configs}
            for rep in range(args.repetitions):
                # Rotate execution order to limit systematic warmup advantages.
                ordered = configs[rep % len(configs):] + configs[:rep % len(configs)]
                for variant, binary, chosen_method in ordered:
                    rows = run(binary, chosen_method, captures, args.timeout_seconds)
                    measures = {
                        "total_ms": sum(float(r["total_ms"]) for r in rows),
                        "full_setup_ms": sum(float(r["full_setup_ms"]) for r in rows),
                        "numeric_update_ms": sum(float(r["numeric_update_ms"]) for r in rows),
                        "krylov_ms": sum(float(r["krylov_ms"]) for r in rows),
                        "iterations_total": sum(int(r["iterations"]) for r in rows),
                        "operator_complexity": float(rows[-1]["amg_operator_complexity"]),
                        "amg_levels": int(rows[-1]["amg_levels"]),
                        "max_true_residual": max(float(r["true_relative_residual"]) for r in rows),
                    }
                    grouped[variant].append(measures)
                    results.append(dict(mesh=mesh, method=method,
                                        variant=variant, repetition=rep, **measures))
            for variant, observations in grouped.items():
                item = dict(mesh=mesh, method=method, variant=variant)
                for field in observations[0]:
                    if field in ("iterations_total", "amg_levels"):
                        item[field] = observations[0][field]
                    elif field == "max_true_residual":
                        item[field] = max(obs[field] for obs in observations)
                    else:
                        item[field] = statistics.median(obs[field] for obs in observations)
                item["reference_total_ms"] = statistics.median(
                    obs["total_ms"] for obs in grouped["v3_map"])
                item["speedup_vs_v3"] = (item["reference_total_ms"] / item["total_ms"]
                                          if item["total_ms"] > 0 else float("nan"))
                item["qualification"] = "PASS_LINEAR_TRUE_RESIDUAL_ONLY"
                summary.append(item)
                print(f"{mesh} {method} {variant}: median={item['total_ms']:.2f}ms "
                      f"speedup={item['speedup_vs_v3']:.2f}x "
                      f"operator_complexity={item['operator_complexity']:.3f} "
                      f"iters={item['iterations_total']}", flush=True)
    for path, rows in ((args.raw, results), (args.summary, summary)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print("NACA_AMG_CSR_QUALIFICATION: PASS_TRUE_RESIDUAL_ONLY", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, ValueError, subprocess.TimeoutExpired) as exc:
        print(f"NACA_AMG_CSR_QUALIFICATION: FAIL {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
