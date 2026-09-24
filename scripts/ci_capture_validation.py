#!/usr/bin/env python3
"""Capture the complete numerical output of the CI validation suite.

This is diagnostic-only: it never changes the CI gate and never stops after
the first failed test.  The authoritative CTest invocation is repeated once
with verbose logging so every model/case output and every iteration history
printed by the executable is preserved in a single artifact.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    log = args.output_dir / "ctest-verbose.log"
    inventory = args.output_dir / "ctest-inventory.log"

    started = time.time()
    listed = subprocess.run(
        ["ctest", "--test-dir", str(args.build_dir), "-N"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    inventory.write_text(listed.stdout, encoding="utf-8")

    # -V preserves the complete stdout/stderr emitted by each test.  This is
    # deliberately a second invocation rather than parsing the shorter
    # --output-on-failure stream, so successful cases are retained as well.
    proc = subprocess.run(
        ["ctest", "--test-dir", str(args.build_dir), "-V"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    elapsed = time.time() - started
    log.write_text(
        f"$ ctest --test-dir {args.build_dir} -V\n"
        f"return_code={proc.returncode}\n"
        f"elapsed_seconds={elapsed:.6f}\n"
        f"--- output ---\n{proc.stdout}",
        encoding="utf-8",
    )

    metadata = {
        "status": "DIAGNOSTIC_ONLY",
        "ctest_return_code": proc.returncode,
        "elapsed_seconds": elapsed,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "commit": os.environ.get("GITHUB_SHA"),
        "workflow": os.environ.get("GITHUB_WORKFLOW"),
        "run_id": os.environ.get("GITHUB_RUN_ID"),
        "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "ref": os.environ.get("GITHUB_REF"),
    }
    (args.output_dir / "diagnostic_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    # Never mask the authoritative CTest result: this evidence step must be
    # allowed to complete after a failure so the artifact remains available.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
