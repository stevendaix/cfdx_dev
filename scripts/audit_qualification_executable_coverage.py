#!/usr/bin/env python3
"""Audit executable coverage of the CFDX engineering qualification matrix.

This is deliberately different from numerical validation: it prevents a
qualification row from being promoted to READY/PASS without a registered
solver-level executable. PLANNED rows are reported as gaps and do not pass as
validated physics.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

EXPECTED = {
    "LAM-COUETTE": "test_couette_quick",
    "LAM-POISEUILLE": "test_poiseuille_quick",
    "INC-GHIA": "test_ghia_cavity_quick",
    "VER-MMS": "test_mms_scalar_diffusion",
    "FORCE-VMFL036": "test_vmfl036_reference_case",
    "EXT-NACA0012": None,
    "INT-CHANNEL": None,
    "SEP-BFS": None,
    "TUR-FLATPLATE": None,
    "TUR-CHANNEL": None,
    "TUR-NACA0012": None,
    "TH-CONV": "test_thermal_vv",
    "CHT-INTERFACE": "test_cht_validation",
    "RAD-P1": "test_radiation_vv",
    "RAD-DOM": "test_radiation_vv",
    "RAD-S2S": "test_s2s_radiation_vv",
}
REQUIRED_EXECUTABLE_STATUSES = {"READY", "RUNNING", "PASS"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("registry", nargs="?", default="docs/validation/CFDX_QUALIFICATION_REGISTRY.json")
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    args = parser.parse_args()

    registry = json.loads(Path(args.registry).read_text(encoding="utf-8"))
    cases = {case["id"]: case for case in registry["cases"]}

    if not (args.build_dir / "CTestTestfile.cmake").exists():
        print("QUALIFICATION_EXECUTABLE_AUDIT: BLOCKED (build directory is not configured)")
        return 2

    listed = subprocess.run(
        ["ctest", "--test-dir", str(args.build_dir), "-N"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout

    failures = []
    planned = []
    for ident, expected in EXPECTED.items():
        case = cases.get(ident)
        if case is None:
            failures.append(f"{ident}: missing from qualification registry")
            continue
        status = case["status"]
        if expected is None:
            planned.append(f"{ident}: {status} -> no solver-level executable yet")
            continue
        if expected not in listed:
            if status in REQUIRED_EXECUTABLE_STATUSES:
                failures.append(f"{ident}: status={status} but CTest executable '{expected}' is not registered")
            else:
                planned.append(f"{ident}: status={status} -> expected executable '{expected}' not registered")
        else:
            print(f"QUALIFICATION_EXECUTABLE id={ident} status={status} test={expected}")

    print("QUALIFICATION_PLANNED_GAPS")
    for item in planned:
        print(f"  - {item}")

    if failures:
        print("QUALIFICATION_EXECUTABLE_AUDIT: FAIL")
        for item in failures:
            print(f"  - {item}")
        return 1

    print("QUALIFICATION_EXECUTABLE_AUDIT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
