"""Paired baseline/candidate on identical captured pressure matrices.

Uses time reported inside the replay (excludes capture loading), alternates
run order to reduce thermal/cache bias, and fails closed on any linear solve
that does not meet the independently recomputed FP64 true-residual gate.
Not a NACA aerodynamic convergence or physical validation.
"""
import argparse
import csv
import io
import math
from pathlib import Path
import statistics
import subprocess
import sys


def run(executable: Path, method: str, manifests: list[Path], timeout: float) -> list[dict]:
    proc = subprocess.run(
        [str(executable.resolve()), method, "reuse", "1", *(str(p.resolve()) for p in manifests)],
        text=True, capture_output=True, check=False, timeout=timeout,
    )
    if proc.returncode:
        raise RuntimeError(f"Replay failed ({executable}, {method}): {proc.stderr[-1000:]}")
    rows = list(csv.DictReader(io.StringIO(proc.stdout)))
    if len(rows) != len(manifests):
        raise RuntimeError(f"Missing replay rows: {executable}: {len(rows)} != {len(manifests)}")
    for row in rows:
        residual = float(row["true_relative_residual"])
        if (row["status"] != "PASS" or not math.isfinite(residual)
                or residual > 1e-10):
            raise RuntimeError(f"Uncertified result: {executable} {method}: {row}")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--candidate", type=Path, required=True)
    ap.add_argument("--captures", type=Path, required=True)
    ap.add_argument("--raw", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--repetitions", type=int, default=5)
    ap.add_argument("--meshes", nargs="+", choices=("N128", "N256"),
                    default=["N128", "N256"])
    ap.add_argument("--methods", nargs="+", choices=("fgmres_sa_fast", "fgmres_boomer_fast"),
                    default=["fgmres_sa_fast", "fgmres_boomer_fast"])
    ap.add_argument("--timeout-seconds", type=float, default=60.0)
    args = ap.parse_args()
    if not 2 <= args.repetitions <= 30:
        ap.error("repetitions must be between 2 and 30")
    records = []
    summaries = []
    for mesh in args.meshes:
        files = sorted(args.captures.glob(f"{mesh}_portable*.manifest"))
        if len(files) != 3:
            raise RuntimeError(f"Expected exactly 3 production {mesh} captures, got {len(files)}")
        for method in args.methods:
            aggregated: dict[str, list[float]] = {"baseline": [], "candidate": []}
            update_times: dict[str, list[float]] = {"baseline": [], "candidate": []}
            for repeat in range(args.repetitions):
                # Reverse the order in alternate pairs to de-bias cold caches.
                variants = [("baseline", args.baseline), ("candidate", args.candidate)]
                if repeat % 2:
                    variants.reverse()
                for label, path in variants:
                    rows = run(path, method, files, args.timeout_seconds)
                    total = 0.0
                    for row in rows:
                        item = dict(row)
                        item.update(version=label, pair=repeat, mesh=mesh)
                        records.append(item)
                        total += float(item["total_ms"])
                        if float(item["numeric_update_ms"]) > 0:
                            update_times[label].append(float(item["numeric_update_ms"]))
                    aggregated[label].append(total)
            old = statistics.median(aggregated["baseline"])
            new = statistics.median(aggregated["candidate"])
            old_update = statistics.median(update_times["baseline"])
            new_update = statistics.median(update_times["candidate"])
            summaries.append(dict(mesh=mesh, method=method,
                median_baseline_ms=old, median_candidate_ms=new,
                total_change_percent=100.0*(new/old - 1.0),
                median_baseline_update_ms=old_update,
                median_candidate_update_ms=new_update,
                update_change_percent=100.0*(new_update/old_update - 1.0),
                repetitions=args.repetitions, snapshots=len(files), qualification="PASS_LINEAR_ONLY"))
            print(f"{mesh:4s} {method:20s} total {old:.2f}->{new:.2f} ms "
                  f"({100*(new/old-1):+.1f}%), updates {old_update:.2f}->{new_update:.2f} ms "
                  f"({100*(new_update/old_update-1):+.1f}%)", flush=True)
    for path, data in ((args.raw, records), (args.summary, summaries)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as stream:
            # Baseline replay binaries may have fewer diagnostic columns
            # than the candidate. Preserve both schemas instead of losing
            # metrics or failing after a lengthy paired campaign.
            columns = list(dict.fromkeys(key for record in data for key in record))
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(data)
    print("PAIRED_NACA_PRESSURE_REPLAY: PASS_TRUE_RESIDUALS physical_qualification=NOT_ATTEMPTED")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as exc:
        print(f"PAIRED_NACA_PRESSURE_REPLAY: FAIL {exc}", file=sys.stderr)
        raise SystemExit(1) from exc