#!/usr/bin/env python3
"""Capture a non-destructive, per-test CFDX validation evidence bundle for CI.

The diagnostic campaign is deliberately independent from the gating CTest step:
every executable is attempted even when another one fails, and every stdout/stderr
stream is preserved.  This prevents an early failure from hiding the numerical
results produced by later models/cases.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path


def discover_ctest_names(build_dir: Path) -> list[str]:
    proc = subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "-N"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    names: list[str] = []
    for line in proc.stdout.splitlines():
        match = re.match(r"\s*Test\s+#\d+:\s+(.+?)\s*$", line)
        if match:
            names.append(match.group(1))
    return names


def run_one(cmd: list[str], log: Path, cwd: Path | None = None) -> int:
    started = time.time()
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    elapsed = time.time() - started
    log.write_text(
        f"$ {' '.join(cmd)}\n"
        f"return_code={proc.returncode}\n"
        f"elapsed_seconds={elapsed:.6f}\n"
        f"--- output ---\n{proc.stdout}",
        encoding="utf-8",
    )
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # CTest is the authoritative inventory: this automatically expands as new
    # validation models/cases are registered in CMake.
    names = discover_ctest_names(args.build_dir)
    records: list[dict[str, object]] = []

    for index, name in enumerate(names, 1):
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or f"test_{index}"
        log = args.output_dir / f"{index:03d}_{safe}.log"
        rc = run_one(
            ["ctest", "--test-dir", str(args.build_dir), "-R", f"^{re.escape(name)}$",
             "--output-on-failure", "--no-tests=error"],
            log,
        )
        records.append({"name": name, "return_code": rc, "log": log.name})

    # Also preserve the raw executable-level output for validation binaries.
    # This catches detailed model diagnostics that a CTest wrapper may truncate.
    executable_names = [
        "test_phase9_acceptance",
        "test_steady_incompressible_solver",
        "test_ghia_cavity",
        "test_analytical_benchmarks",
        "test_benchmark_matrix",
        "test_level_b_reference_benchmarks",
        "test_level_c_coupled_verification",
        "test_numerical_model_verification",
        "test_m1_m4_validation",
        "test_m2_m4_solver_validation",
        "test_cht_validation",
        "test_fluent_vmfl_reference",
    ]
    exe_dir = args.output_dir / "executables"
    exe_dir.mkdir(exist_ok=True)

    for name in executable_names:
        exe = args.build_dir / name
        if not exe.is_file():
            nested = args.build_dir / "tests" / name
            exe = nested if nested.is_file() else exe
        if not exe.is_file():
            records.append({"name": name, "kind": "executable", "return_code": -1,
                            "status": "MISSING"})
            continue
        log = exe_dir / f"{name}.log"
        rc = run_one([str(exe)], log)
        records.append({"name": name, "kind": "executable", "return_code": rc,
                        "status": "PASS" if rc == 0 else "FAIL", "log": str(log.relative_to(args.output_dir))})

    (args.output_dir / "diagnostic_results.json").write_text(
        json.dumps(
            {
                "build_dir": str(args.build_dir),
                "test_count": len(names),
                "records": records,
                "failed": [r for r in records if r.get("return_code") not in (0, None)],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    # Diagnostic collection must not replace the real CI gate.  Always return
    # success so the following upload step can execute; the preceding/following
    # gating tests remain authoritative.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
