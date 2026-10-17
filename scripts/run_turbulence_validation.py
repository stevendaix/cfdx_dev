#!/usr/bin/env python3
"""Run the executable part of CFDX turbulence qualification without overstating evidence.

This driver runs existing CTest targets and records missing gates as BLOCKED.
It does not create a turbulent physical-validation PASS from unit tests.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REQUIRED_TESTS = (
    "test_turbulence_qualification_equations",
    "test_phase10_turbulence_hardening",
    "test_n11_sa_negative_nu_tilde_policy",
    "test_n11_bound_reporting_hardening",
)
BASELINE_TESTS = ("test_channel_qualification",)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    print("$ " + " ".join(command), flush=True)
    return subprocess.run(command, text=True, check=False)


def registered_tests(build_dir: Path) -> set[str]:
    result = run(["ctest", "--test-dir", str(build_dir), "-N"])
    if result.returncode != 0:
        raise RuntimeError(f"ctest -N failed with exit code {result.returncode}")
    return set(re.findall(r"Test\s+#\d+:\s+(\S+)", result.stdout))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/turbulence-qualification/report.json"),
    )
    args = parser.parse_args()
    build_dir = args.build_dir.resolve()
    if not (build_dir / "CTestTestfile.cmake").is_file():
        print(f"ERROR: build directory is not configured: {build_dir}", file=sys.stderr)
        return 2

    try:
        available = registered_tests(build_dir)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    records: list[dict[str, object]] = []
    missing_required: list[str] = []
    failed = False
    for category, names in (("turbulence", REQUIRED_TESTS), ("laminar_baseline", BASELINE_TESTS)):
        for name in names:
            if name not in available:
                status = "BLOCKED_NOT_REGISTERED"
                if category == "turbulence":
                    missing_required.append(name)
                records.append({"category": category, "test": name, "status": status})
                print(f"{status}: {name}")
                continue
            result = run([
                "ctest", "--test-dir", str(build_dir), "-R", f"^{re.escape(name)}$",
                "--output-on-failure",
            ])
            status = "PASS" if result.returncode == 0 else "FAIL"
            failed |= result.returncode != 0
            records.append({
                "category": category,
                "test": name,
                "status": status,
                "exit_code": result.returncode,
            })
            print(f"{status}: {name}")

    report = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "build_dir": str(build_dir),
        "scope": "existing component/equation tests plus laminar baseline; not physical validation",
        "results": records,
        "physical_validation": {
            "status": "BLOCKED_NOT_IMPLEMENTED_AS_QUALIFIED_CAMPAIGN",
            "required_cases": ["TUR-CHANNEL", "TUR-FLATPLATE", "TUR-NACA0012"],
            "note": "A unit-test pass must not be promoted to turbulent solver or physical-validation PASS.",
        },
        "turbulence_transport_mms": {
            "status": "BLOCKED_NOT_REGISTERED_AS_DEDICATED_CAMPAIGN",
            "required_models": ["k-epsilon", "k-omega SST", "Spalart-Allmaras"],
            "note": "The Phase 13 acceptance plan identifies isolated transport MMS as missing evidence.",
        },
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Report: {output}")
    if failed:
        return 1
    if missing_required:
        return 2
    # Exit 2 is intentional while physical/MMS gates remain blocked.
    print("BLOCKED: component tests completed, but turbulence qualification is not closed.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
