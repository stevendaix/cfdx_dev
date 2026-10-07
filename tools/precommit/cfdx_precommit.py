"""Fast, change-aware CFDX pre-commit gate.

This is not numerical qualification. It checks that changed additions can still
configure/build, that cheap relevant tests execute, and that obvious test
bypasses were not introduced. Long validation and N1-N17 remain CI gates.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("CFDX_PRECOMMIT_BUILD_DIR", ROOT / "build-precommit"))
JOBS = os.environ.get("CFDX_PRECOMMIT_JOBS", "4")

SMOKE_TESTS = [
    "test_convergence_history",
    "test_execution_summary",
    "test_hdf5_roundtrip",
    "test_case_hdf5_io",
]
CPP_EXT = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx"}


def run(cmd: list[str], *, cwd: Path = ROOT) -> None:
    print("+", " ".join(cmd), flush=True)
    result = subprocess.run(cmd, cwd=cwd, check=False)
    if result.returncode:
        raise SystemExit(result.returncode)


def output(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, cwd=ROOT, text=True, stderr=subprocess.STDOUT)


def changed_files() -> list[str]:
    range_spec = os.environ.get("CFDX_PRECOMMIT_RANGE")
    if range_spec:
        return [
            x for x in output(
                ["git", "diff", range_spec, "--name-only", "--diff-filter=ACMR"]
            ).splitlines()
            if x
        ]
    return [
        x for x in output(
            ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"]
        ).splitlines()
        if x
    ]


def changed_diff() -> str:
    range_spec = os.environ.get("CFDX_PRECOMMIT_RANGE")
    if range_spec:
        return output(["git", "diff", range_spec, "--unified=0", "--"])
    return output(["git", "diff", "--cached", "--unified=0", "--"])


def existing_ctest_names() -> set[str]:
    if not (BUILD / "CTestTestfile.cmake").exists():
        return set()
    listing = output(["ctest", "--test-dir", str(BUILD), "-N"])
    return set(re.findall(r"Test #\d+: ([^\n]+)", listing))


def configure() -> None:
    BUILD.mkdir(parents=True, exist_ok=True)
    if (BUILD / "CMakeCache.txt").exists():
        return
    cmd = [
        "cmake", "-S", ".", "-B", str(BUILD),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCFDX_BUILD_TESTS=ON",
        "-DCFDX_BUILD_RUNTIME_TESTS=ON",
        "-DCFDX_ENABLE_GPU=OFF",
        "-DCFDX_ENABLE_LONG_VALIDATION=OFF",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
    ]
    if shutil.which("ninja"):
        cmd[6:6] = ["-G", "Ninja"]
    run(cmd)


def build_targets(targets: list[str]) -> None:
    if targets:
        run(["cmake", "--build", str(BUILD), "--parallel", JOBS, "--target", *targets])


def run_ctest(names: list[str]) -> None:
    available = existing_ctest_names()
    selected = [name for name in names if name in available]
    if not selected:
        print("CFDX pre-commit: no smoke test registered; build-only guard remains active.")
        return
    for name in selected:
        run([
            "ctest", "--test-dir", str(BUILD),
            "-R", f"^{re.escape(name)}$", "--output-on-failure",
        ])


def python_guard(files: list[str]) -> None:
    py = [x for x in files if x.endswith(".py")]
    if not py:
        return
    run([sys.executable, "-m", "compileall", "-q", *py])
    if shutil.which("ruff"):
        run(["ruff", "check", *py])
    else:
        print("CFDX pre-commit: ruff not installed; Python compile check still ran.")


def integrity_guard(diff: str) -> None:
    added: list[str] = []
    current_file = ""
    for line in diff.splitlines():
        if line.startswith("+++ b/"):
            current_file = line[6:]
            continue
        if (
            line.startswith("+")
            and not line.startswith("+++")
            and current_file != "tools/precommit/cfdx_precommit.py"
        ):
            added.append(line[1:])
    forbidden = [
        (r"\bDISABLED_[A-Za-z0-9_]+", "new disabled test"),
        (r"\bGTEST_SKIP\s*\(", "new GTEST_SKIP"),
        (r"\bSKIP_RETURN_CODE\b", "new SKIP_RETURN_CODE"),
    ]
    hits = [
        label for pattern, label in forbidden
        if any(re.search(pattern, line) for line in added)
    ]
    if hits:
        print("CFDX pre-commit: test-integrity guard rejected staged changes:", file=sys.stderr)
        for hit in sorted(set(hits)):
            print(f"  - {hit}", file=sys.stderr)
        print(
            "Make intentional skips explicit in the PR/CI review; do not bypass "
            "a test merely to obtain a green local commit.",
            file=sys.stderr,
        )
        raise SystemExit(1)


def main() -> int:
    files = changed_files()
    if not files:
        return 0

    print(f"CFDX pre-commit: {len(files)} changed files")
    integrity_guard(changed_diff())

    python_changed = any(x.endswith(".py") for x in files)
    cpp_changed = any(Path(x).suffix.lower() in CPP_EXT for x in files)
    cmake_changed = any(
        x == "CMakeLists.txt" or x.startswith("cmake/") or x.endswith(".cmake")
        for x in files
    )

    if python_changed:
        python_guard(files)

    if cpp_changed or cmake_changed:
        configure()
        available_smoke = existing_ctest_names()
        targets = ["cfdx_core", *(name for name in SMOKE_TESTS if name in available_smoke)]
        for path in files:
            parts = Path(path).parts
            if len(parts) >= 3 and parts[0] == "tests" and Path(path).suffix.lower() in CPP_EXT:
                stem = Path(path).stem
                if stem.startswith("test_"):
                    targets.append(stem)
        build_targets(list(dict.fromkeys(targets)))
        run_ctest(SMOKE_TESTS)

    print("CFDX pre-commit: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
