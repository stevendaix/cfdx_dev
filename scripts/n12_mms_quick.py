#!/usr/bin/env python3
"""Run the currently implemented N12 MMS quick campaign and retain evidence.

This is an orchestration tool, not a numerical oracle.  It runs existing
authoritative CTest MMS executables, records their exit status, timing and
output, and fails if an expected test is missing or fails.  It never changes
solver tolerances and never converts a timeout/failure into PASS.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

TESTS = (
    "test_mms_scalar_diffusion",
    "test_transient_mms",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    inventory = subprocess.run(
        ["ctest", "--test-dir", str(args.build_dir), "-N"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    inventory_text = inventory.stdout
    (args.output_dir / "ctest-inventory.log").write_text(
        inventory_text, encoding="utf-8"
    )

    missing = [name for name in TESTS if name not in inventory_text]
    results: list[dict[str, object]] = []

    for name in TESTS:
        record: dict[str, object] = {"name": name}
        if name in missing:
            record.update({"status": "MISSING", "return_code": None})
            results.append(record)
            continue

        log_path = args.output_dir / f"{name}.log"
        started = time.monotonic()
        try:
            proc = subprocess.run(
                [
                    "ctest",
                    "--test-dir",
                    str(args.build_dir),
                    "--tests-regex",
                    f"^{name}$",
                    "--no-tests=error",
                    "--output-on-failure",
                    "-V",
                    "--timeout",
                    str(args.timeout),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            elapsed = time.monotonic() - started
            log_path.write_text(proc.stdout, encoding="utf-8")
            status = "PASS" if proc.returncode == 0 else "FAIL"
            record.update(
                {
                    "status": status,
                    "return_code": proc.returncode,
                    "elapsed_seconds": round(elapsed, 3),
                    "log": log_path.name,
                }
            )
        except subprocess.TimeoutExpired as exc:
            elapsed = time.monotonic() - started
            output = exc.stdout or ""
            log_path.write_text(output, encoding="utf-8")
            record.update(
                {
                    "status": "TIMEOUT",
                    "return_code": None,
                    "elapsed_seconds": round(elapsed, 3),
                    "log": log_path.name,
                }
            )
        results.append(record)

    status = "PASS" if inventory.returncode == 0 and all(
        item["status"] == "PASS" for item in results
    ) else "FAIL"

    report = {
        "status": status,
        "campaign": "N12-MMS-QUICK",
        "tests": list(TESTS),
        "timeout_seconds": args.timeout,
        "results": results,
    }
    (args.output_dir / "n12-mms-quick.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    print(json.dumps(report, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
