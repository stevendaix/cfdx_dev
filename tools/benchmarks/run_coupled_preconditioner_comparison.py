"""CFDX V5 fair production-captured COUPLED benchmarking, fail-closed.

Each grid is captured ONCE at the third nonlinear iteration of the real
Ghia Re=100 solver. Five independent processes replay the exact snapshot.
A different five-process campaign measures full physical convergence.

This driver does not infer convergence from Krylov status: it checks the
independently computed true residual and the physical solver's final gates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

METHODS = ("amg", "schwarz", "bfbt", "lsc", "pcd")
KEY_VALUE = re.compile(r"([a-z][a-z0-9_]*)=([^\s]+)")
EXPECTED = {"linear": "CFDX_LINEAR_V5 ", "physical": "CFDX_PHYSICAL_V5 "}


def run_one(binary: Path, args: list[str], timeout: float) -> dict:
    start = time.perf_counter()
    try:
        proc = subprocess.run(
            [str(binary), *args],
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        raw = proc.stdout + proc.stderr
        duration = time.perf_counter() - start
        return {
            "exit_code": proc.returncode,
            "process_elapsed_s": duration,
            "raw": raw,
            "timed_out": False,
        }
    except subprocess.TimeoutExpired as exc:
        partial = exc.stdout or b""
        if isinstance(partial, bytes):
            partial = partial.decode("utf8", "replace")
        return {
            "exit_code": 124,
            "process_elapsed_s": time.perf_counter() - start,
            "raw": partial + "\nTIMEOUT",
            "timed_out": True,
        }


def parse_result(run: dict, kind: str) -> dict:
    prefix = EXPECTED[kind]
    lines = [line for line in run["raw"].splitlines() if line.startswith(prefix)]
    fields = dict(KEY_VALUE.findall(lines[-1])) if len(lines) == 1 else {}
    return {
        **fields,
        "exit_code": run["exit_code"],
        "process_elapsed_s": round(run["process_elapsed_s"], 6),
        "timed_out": run["timed_out"],
    }


def quality(result: dict, kind: str, expected_hash: str = "") -> tuple[bool, str]:
    if result.get("exit_code") != 0:
        return False, "execution_failed_or_timed_out"
    if result.get("converged") not in ("1", "true", "True"):
        return False, "solver_not_converged"
    if kind == "linear":
        if result.get("fingerprint") != expected_hash:
            return False, "capture_fingerprint_mismatch"
        try:
            rr = float(result["true_relative_residual"])
        except (KeyError, ValueError):
            return False, "missing_true_residual"
        if not 0 <= rr <= 1e-10:
            return False, "linear_true_residual_gate_failed"
    else:
        try:
            momentum = float(result["momentum_relative"])
            continuity = float(result["continuity_normalized"])
        except (KeyError, ValueError):
            return False, "missing_physical_metrics"
        if not 0 <= momentum <= 1e-8 or not 0 <= continuity <= 1e-8:
            return False, "physical_convergence_gate_failed"
    return True, "passed"


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_table(path: Path, records: list[dict]) -> None:
    keys = sorted(set().union(*(r.keys() for r in records)))
    with path.open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(records)


def test_protocol() -> None:
    sample = {
        "exit_code": 0,
        "converged": "1",
        "fingerprint": "abc",
        "true_relative_residual": "8e-11",
    }
    assert quality(sample, "linear", "abc")[0]
    assert not quality(sample, "linear", "xyz")[0]
    assert not quality(
        {**sample, "true_relative_residual": "1.2e-10"}, "linear", "abc"
    )[0]
    assert not quality({**sample, "exit_code": 124}, "linear", "abc")[0]
    assert not quality({**sample, "converged": "0"}, "linear", "abc")[0]
    assert not quality({"exit_code": 0, "converged": "1"}, "physical")[0]
    phys = {
        "exit_code": 0,
        "converged": "1",
        "momentum_relative": "1e-9",
        "continuity_normalized": "1e-9",
    }
    assert quality(phys, "physical")[0]
    assert not quality({**phys, "continuity_normalized": "2e-8"}, "physical")[0]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path)
    p.add_argument("--output", type=Path, default=Path("coupled_comparison_v5"))
    p.add_argument("--grids", nargs="+", type=int, default=[6, 12])
    p.add_argument("--outer-budget", type=int, default=350)
    p.add_argument("--timeout", type=float, default=120)
    p.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="independent processes for reproducible median timings",
    )
    p.add_argument("--skip-physical", action="store_true")
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args()
    test_protocol()
    if args.selftest:
        print("protocol selftest passed: 8 fail-closed checks")
        return 0
    if args.binary is None or not args.binary.is_file():
        p.error("provide a compiled --binary")
    if any(n < 4 or n > 64 for n in args.grids):
        p.error("supported grids: 4..64")
    if (
        args.outer_budget < 1
        or args.timeout <= 0
        or args.repeats < 1
        or args.repeats > 20
    ):
        p.error("invalid budget")
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = {
        "binary_sha256": hash_file(args.binary),
        "binary": str(args.binary.resolve()),
        "system": platform.platform(),
        "python": sys.version,
        "cpu_count": os.cpu_count(),
        "grids": args.grids,
        "outer_budget": args.outer_budget,
        "repeats": args.repeats,
        "linear_true_residual_gate": 1e-10,
        "nonlinear_relative_and_continuity_gates": 1e-8,
        "snapshot_iteration": 3,
        "measurement": (
            "isolated process per variant; compare VmRSS and glibc live-heap deltas. "
            "ru_maxrss is NOT compared because forked Python processes inherit "
            "their parent's historic high-water mark"
        ),
    }
    records = []
    invalid = []
    for grid in args.grids:
        snapshot = args.output / f"ghia_re100_{grid}x{grid}_stage3.cfdxbin"
        capture = run_one(
            args.binary, ["capture", str(grid), str(snapshot)], args.timeout
        )
        (args.output / f"capture_{grid}.log").write_text(
            capture["raw"], encoding="utf8"
        )
        lines = [
            line
            for line in capture["raw"].splitlines()
            if line.startswith("CFDX_CAPTURE_V5 ")
        ]
        captured = dict(KEY_VALUE.findall(lines[-1])) if len(lines) == 1 else {}
        fingerprint = captured.get("fingerprint", "")
        if capture["exit_code"] != 0 or not fingerprint or not snapshot.is_file():
            invalid.append(f"capture {grid}x{grid} failed")
            print(f"ERROR: capture {grid}x{grid} failed", file=sys.stderr)
            continue
        print(
            f"CAPTURE grid={grid} fingerprint={fingerprint} sha256={hash_file(snapshot)[:16]}",
            flush=True,
        )
        for method in METHODS:
            for kind in ("linear", "physical"):
                if kind == "physical" and args.skip_physical:
                    continue
                argv = (
                    ["replay", method, str(snapshot)]
                    if kind == "linear"
                    else ["physical", method, str(grid), str(args.outer_budget)]
                )
                for rep in range(1, args.repeats + 1):
                    executed = run_one(args.binary, argv, args.timeout)
                    (args.output / f"{kind}_{grid}_{method}_r{rep}.log").write_text(
                        executed["raw"], encoding="utf8"
                    )
                    row = parse_result(executed, kind)
                    row.update(kind=kind, grid=grid, method=method, repeat=rep)
                    # Linux inherited ru_maxrss can reflect Python's pre-exec
                    # memory rather than CFDX's peak. It is not a valid cross-
                    # method memory comparison; preserve it only in raw logs.
                    row.pop("process_peak_rss_kib", None)
                    passed, reason = quality(row, kind, fingerprint)
                    row["qualified"] = passed
                    row["gate_reason"] = reason
                    records.append(row)
                    if not passed:
                        invalid.append(
                            f"{kind} {grid}x{grid} {method} r{rep}: {reason}"
                        )
                    print(
                        f"{kind.upper()} grid={grid} method={method} rep={rep} "
                        f"pass={int(passed)} "
                        f"iters={row.get('outer_iterations', row.get('iterations', '?'))} "
                        f"true_rr={row.get('true_relative_residual', '?')} "
                        f"wall_ms={row.get('total_linear_ms', row.get('wall_ms', '?'))}",
                        flush=True,
                    )
    # Cross-method physical consistency is independent of convergence flags:
    # compare finished velocity norms and the pressure range on each grid.
    if not args.skip_physical:
        for grid in args.grids:
            rows = [r for r in records if r["grid"] == grid and r["kind"] == "physical"]
            ref = next(
                (r for r in rows if r["method"] == "amg" and r["qualified"]), None
            )
            if ref is None or len(rows) != len(METHODS) * args.repeats:
                invalid.append(
                    f"physical {grid}x{grid}: incomplete physical comparison"
                )
                continue
            for item in rows:
                try:
                    for key in ("velocity_rms_u", "velocity_rms_v", "pressure_range"):
                        if abs(float(item[key]) - float(ref[key])) > 1e-7:
                            raise ValueError(
                                f"{key} differs from common final solution"
                            )
                    if abs(float(item["boundary_mass_flux"])) > 1e-8:
                        raise ValueError("global boundary flux not conserved")
                except (KeyError, ValueError) as exc:
                    item["qualified"] = False
                    item["gate_reason"] = f"physical_state_mismatch: {exc}"
                    invalid.append(f"physical {grid}x{grid} {item['method']}: {exc}")
    write_table(args.output / "results.csv", records)
    (args.output / "results.json").write_text(
        json.dumps(
            {
                "metadata": metadata,
                "results": records,
                "invalid": invalid,
                "qualified": not invalid,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf8",
    )
    print(
        "FINAL",
        "qualified" if not invalid else "FAILED",
        "invalid_count=",
        len(invalid),
    )
    return 0 if not invalid else 2


if __name__ == "__main__":
    raise SystemExit(main())
