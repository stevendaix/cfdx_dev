#!/usr/bin/env python3
"""Run the complete M1-M4 verification campaign (Levels A, B and C)."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


LEVELS = {
    "A": ["test_analytical_benchmarks", "test_benchmark_matrix"],
    "B": ["test_level_b_reference_benchmarks"],
    "C": ["test_level_c_coupled_verification"],
}


def find_executable(build_dir: Path, name: str) -> Path | None:
    for p in (build_dir / name, build_dir / "tests" / name):
        if p.is_file():
            return p
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    args = parser.parse_args()

    failed = False
    print("CFDX M1-M4 VERIFICATION")
    print("=======================")

    for level, executables in LEVELS.items():
        print(f"\nLEVEL {level}")
        for name in executables:
            exe = find_executable(args.build_dir, name)
            if exe is None:
                print(f"  {name}: NOT BUILT", file=sys.stderr)
                failed = True
                continue
            result = subprocess.run([str(exe)], check=False)
            status = "PASS" if result.returncode == 0 else "FAIL"
            print(f"  {name}: {status}")
            failed |= result.returncode != 0

    print("\nOVERALL")
    print("=======")
    print("STATUS: FAIL" if failed else "STATUS: PASS")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
