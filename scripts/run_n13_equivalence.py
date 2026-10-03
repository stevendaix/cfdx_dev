#!/usr/bin/env python3
"""Run the N13 cross-backend numerical-equivalence campaign."""
from __future__ import annotations
import argparse
import json
import subprocess
from pathlib import Path

CAMPAIGNS = {
    "serial_mpi": ["test_n13_mpi_equivalence"],
    "assembled_matrix_free": ["test_matrix_free_fv_operator", "test_backend_equivalence"],
    "deterministic_reduction": ["test_n13_mpi_equivalence", "test_mpi_deterministic"],
    "restart": ["test_n13_restart_equivalence", "test_restart_mapping"],
}

def run_test(build_dir: Path, name: str) -> tuple[str, int]:
    proc = subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "-R", f"^{name}$", "--output-on-failure"],
        check=False,
    )
    return name, proc.returncode

def run_cuda_test(build_dir: Path) -> tuple[str, int]:
    exe = build_dir / "test_n13_cuda_equivalence"
    if not exe.exists():
        return "test_n13_cuda_equivalence", 2
    proc = subprocess.run([str(exe)], check=False)
    return "test_n13_cuda_equivalence", proc.returncode

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--json", type=Path, default=Path("artifacts/n13_equivalence.json"))
    args = parser.parse_args()
    if not (args.build_dir / "CTestTestfile.cmake").exists():
        raise SystemExit(f"build directory is not configured: {args.build_dir}")

    report = {
        "schema_version": "1.0",
        "package": "N13",
        "name": "Cross-backend numerical equivalence",
        "tolerances": {
            "operator_l2_relative": 1e-12,
            "operator_linf_relative": 1e-12,
            "restart": 0.0,
            "deterministic_reduction": 1e-14,
        },
        "campaigns": {},
    }
    failures = 0
    for campaign, tests in CAMPAIGNS.items():
        entries = []
        for test in tests:
            name, rc = run_test(args.build_dir, test)
            status = "PASS" if rc == 0 else "FAIL"
            entries.append({"test": name, "status": status, "exit_code": rc})
            if status == "FAIL":
                failures += 1
        report["campaigns"][campaign] = entries

    name, rc = run_cuda_test(args.build_dir)
    report["campaigns"]["cuda"] = [{
        "test": name,
        "status": "PASS" if rc == 0 else ("BLOCKED" if rc == 2 else "FAIL"),
        "exit_code": rc,
    }]
    if rc not in (0, 2):
        failures += 1

    report["closure"] = "PASS" if failures == 0 else "FAIL"
    report["note"] = "CUDA is a separate hardware gate; unavailable CUDA hardware is BLOCKED, never PASS."
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
