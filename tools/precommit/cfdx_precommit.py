#!/usr/bin/env python3
"""Fast, change-aware CFDX pre-commit gate.

This is deliberately not a numerical qualification suite. It answers:
1. Can the changed C++/CMake code still configure and compile?
2. Can the most relevant cheap tests execute?
3. Do changed Python modules compile/lint and do focused tests pass?
4. Did the staged diff introduce obvious test-integrity bypasses?

Long validation, CUDA qualification, MPI campaigns, N1-N17 closure and
validation-total remain CI responsibilities.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("CFDX_PRECOMMIT_BUILD_DIR", ROOT / "build-precommit"))
JOBS = os.environ.get("CFDX_PRECOMMIT_JOBS", "4")

SMOKE_TESTS = [
    "test_validation_simple",
    "test_convergence_history",
    "test_execution_summary",
    "test_hdf5_roundtrip",
    "test_case_hdf5_io",
]

CPP_EXT = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx"}

def run(cmd: list[str], *, cwd: Path = ROOT) -> None:
    print("+", " ".join(cmd))
    p = subprocess.run(cmd, cwd=cwd)
    if p.returncode:
        raise SystemExit(p.returncode)

def output(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, cwd=ROOT, text=True, stderr=subprocess.STDOUT)

def staged_files() -> list[str]:
    return [x for x in output(["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"]).splitlines() if x]

def staged_diff() -> str:
    return output(["git", "diff", "--cached", "--unified=0", "--"])

def existing_ctest_names() -> set[str]:
    if not (BUILD / "CTestTestfile.cmake").exists():
        return set()
    text = output(["ctest", "--test-dir", str(BUILD), "-N"])
    return set(re.findall(r"Test #\d+: ([^\\n]+)", text))

def configure() -> None:
    BUILD.mkdir(parents=True, exist_ok=True)
    if (BUILD / "CMakeCache.txt").exists():
        return
    run([
        "cmake", "-S", ".", "-B", str(BUILD), "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCFDX_BUILD_TESTS=ON",
        "-DCFDX_BUILD_RUNTIME_TESTS=ON",
        "-DCFDX_ENABLE_GPU=OFF",
        "-DCFDX_ENABLE_LONG_VALIDATION=OFF",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
    ])

def build_targets(targets: list[str]) -> None:
    if not targets:
        return
    run(["cmake", "--build", str(BUILD), "--parallel", JOBS, "--target", *targets])

def run_ctest(names: list[str]) -> None:
    available = existing_ctest_names()
    selected = [x for x in names if x in available]
    if not selected:
        print("CFDX pre-commit: no selected smoke tests are registered; build-only guard remains active.")
        return
    for name in selected:
        run(["ctest", "--test-dir", str(BUILD), "-R", f"^{re.escape(name)}$", "--output-on-failure"])

def python_guard(files: list[str]) -> None:
    py = [x for x in files if x.endswith(".py")]
    if not py:
        return
    run([sys.executable, "-m", "compileall", "-q", *py])
    try:
        run([sys.executable, "-m", "ruff", "check", *py])
    except SystemExit as exc:
        if exc.code != 0:
            print("CFDX pre-commit: ruff is required for changed Python files.", file=sys.stderr)
            raise

def integrity_guard(diff: str) -> None:
    added = [line[1:] for line in diff.splitlines() if line.startswith("+") and not line.startswith("+++")]
    forbidden = [
        (r"\\bDISABLED_[A-Za-z0-9_]+", "new disabled test"),
        (r"\\bGTEST_SKIP\\s*\\(", "new unconditional GTEST_SKIP"),
        (r"\\bSKIP_RETURN_CODE\\b", "new SKIP_RETURN_CODE"),
    ]
    hits = []
    for pattern, label in forbidden:
        if any(re.search(pattern, line) for line in added):
            hits.append(label)
    if hits:
        print("CFDX pre-commit: test-integrity guard rejected staged changes:", file=sys.stderr)
        for hit in sorted(set(hits)):
            print(f"  - {hit}", file=sys.stderr)
        print("If a skip is genuinely required, make the exception explicit in the PR/CI review rather than bypassing the guard.", file=sys.stderr)
        raise SystemExit(1)

def main() -> int:
    files = staged_files()
    if not files:
        return 0

    print(f"CFDX pre-commit: {len(files)} staged files")
    integrity_guard(staged_diff())

    py = [x for x in files if x.endswith(".py")]
    cpp = [x for x in files if Path(x).suffix.lower() in CPP_EXT]
    cmake_changed = any(
        x == "CMakeLists.txt" or x.startswith("cmake/") or x.endswith(".cmake")
        for x in files
    )

    if py:
        python_guard(files)

    if cpp or cmake_changed:
        configure()
        targets = ["cfdx_core"]
        for path in files:
            p = Path(path)
            if len(p.parts) >= 3 and p.parts[0] == "tests" and p.suffix.lower() in CPP_EXT:
                stem = p.stem
                if stem.startswith("test_"):
                    targets.append(stem)
        build_targets(list(dict.fromkeys(targets)))
        run_ctest(SMOKE_TESTS)

    print("CFDX pre-commit: PASS")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
