#!/usr/bin/env python3
"""Run the deterministic M1-M4 analytical verification executables.

Usage:
    python scripts/run_validation.py --build-dir build
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    args = parser.parse_args()

    candidates = [
        args.build_dir / "test_analytical_benchmarks",
        args.build_dir / "tests" / "test_analytical_benchmarks",
    ]
    executable = next((p for p in candidates if p.is_file()), None)
    if executable is None:
        print(f"Validation executable not found under {args.build_dir}", file=sys.stderr)
        print("Configure and build CFDX first.", file=sys.stderr)
        return 2

    completed = subprocess.run([str(executable)], check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
