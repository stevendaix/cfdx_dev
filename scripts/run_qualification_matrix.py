"""Explicit, fail-closed CFDX CTest qualification runner.

This tool runs the actual CFDX binaries and does not reimplement physics.
A CTest PASS means a test's own gates passed; model-level qualification
requires all associated independent physical gates and refinements.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass(frozen=True)
class Case:
    name: str
    target: str
    level: str  # component, numerical, physical_prequalification, physical
    timeout: int


QUICK = (
    Case(
        "test_thermal_forced_convection_channel",
        "test_thermal_forced_convection_channel",
        "numerical",
        120,
    ),
    Case(
        "test_natural_convection_de_vahl_davis",
        "test_natural_convection_de_vahl_davis",
        "physical_prequalification",
        450,
    ),
    Case("test_thermal_vv", "test_thermal_vv", "component", 90),
    Case("test_radiation_vv", "test_radiation_vv", "component", 90),
    Case("test_s2s_radiation_vv", "test_s2s_radiation_vv", "component", 90),
    Case(
        "test_production_s2s_energy_outer_loop_vv",
        "test_production_s2s_energy_outer_loop_vv",
        "component",
        90,
    ),
    Case(
        "test_thermal_radiation_model_matrix",
        "test_thermal_radiation_model_matrix",
        "component",
        90,
    ),
    Case(
        "test_turbulence_qualification_equations",
        "test_turbulence_qualification_equations",
        "component",
        90,
    ),
    Case(
        "test_production_turbulence_contract",
        "test_production_turbulence_contract",
        "component",
        90,
    ),
)

EXTENDED = (
    Case("test_ghia_cavity_quick", "test_ghia_cavity", "numerical", 180),
    Case(
        "test_n9_s9_rhie_chow_consistency",
        "test_n9_s9_rhie_chow_consistency",
        "numerical",
        180,
    ),
)

# Intentionally separate: full campaigns must not be silently discarded
# when only the quick subset passes. Long tests must be registered at CMake
# configuration time (-DCFDX_ENABLE_LONG_VALIDATION=ON).
LONG = (
    Case(
        "test_naca0012_qualification",
        "test_naca0012_qualification",
        "physical",
        10800,
    ),
    Case(
        "test_naca0012_wall_shear_audit",
        "test_naca0012_wall_shear_audit",
        "physical",
        10800,
    ),
    Case("test_ghia_cavity_profile_study", "test_ghia_cavity", "physical", 5400),
)


def check_ctest_output(returncode: int, output: str, test_name: str) -> bool:
    """Return success only if CTest actually executed exactly one test."""
    # Do not mistake `No tests were found` or absent binaries for a success.
    summary = re.search(
        r"100% tests passed, 0 tests failed out of (\d+)\b", output
    )
    # CTest can auto-add FIXTURES_SETUP / FIXTURES_CLEANUP cases, so the
    # requested test need not be the only one in a successful run.
    target = re.search(
        r"\b\d+/\d+ Test\s+#\d+:\s+"
        + re.escape(test_name)
        + r"\s+\.{2,}\s+Passed\b",
        output,
    )
    return bool(
        returncode == 0
        and summary
        and int(summary.group(1)) >= 1
        and target
        and "No tests were found" not in output
    )


def get_registered_tests(build_dir: Path) -> set[str]:
    p = subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "--show-only=json-v1"],
        capture_output=True,\n        text=True,\n        check=True,\n        timeout=60,
    )
    return {entry["name"] for entry in json.loads(p.stdout)["tests"]}


def run_case(case: Case, build_dir: Path, logs: Path, jobs: int) -> dict:
    row = asdict(case)
    row.update(status="NOT_RUN", reason="", elapsed_s=None)
    started = time.monotonic()
    try:
        build = subprocess.run(
            [
                "cmake",
                "--build",
                str(build_dir),
                "--target",
                case.target,
                "-j",
                str(jobs),
            ],
            capture_output=True,\n            text=True,\n            timeout=600,\n            check=False,
        )
        (logs / f"{case.name}.build.log").write_text(\n            build.stdout + build.stderr\n        )
        if build.returncode != 0:
            row.update(status="BUILD_FAIL", reason=f"build returned {build.returncode}")
            return row
        completed = subprocess.run(
            [
                "ctest",
                "--test-dir",
                str(build_dir),
                "-R",
                f"^{re.escape(case.name)}$",
                "-V",
                "--timeout",
                str(case.timeout),
            ],
            capture_output=True,\n            text=True,\n            timeout=case.timeout + 90,\n            check=False,
        )
        output = completed.stdout + completed.stderr
        (logs / f"{case.name}.ctest.log").write_text(output)
        row["status"] = (
            "PASS"
            if check_ctest_output(completed.returncode, output, case.name)
            else "FAIL"
        )
        row["reason"] = (
            ""
            if row["status"] == "PASS"
            else f"ctest exit {completed.returncode}; see log"
        )
    except subprocess.TimeoutExpired as exc:
        row.update(status="TIMEOUT", reason=f"subprocess timeout: {exc}")
    except OSError as exc:
        row.update(status="INFRA_FAIL", reason=str(exc))
    finally:
        row["elapsed_s"] = round(time.monotonic() - started, 3)
    return row


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument(\n        "--tier", choices=("quick", "extended", "full"), default="quick"\n    )
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument(
        "--preflight-only",
        action="store_true",
        help="Only check CTest registration; do not run or build",
    )
    args = p.parse_args(argv)
    if args.jobs < 1 or args.jobs > 16:
        p.error("--jobs must be between 1 and 16")
    build_dir = args.build_dir.resolve()
    if not (build_dir / "CMakeCache.txt").is_file():
        p.error("--build-dir must contain an existing CMakeCache.txt")
    logs = args.output_dir.resolve()
    logs.mkdir(parents=True, exist_ok=True)
    cases = (
        QUICK
        + (EXTENDED if args.tier in ("extended", "full") else ())
        + (LONG if args.tier == "full" else ())
    )
    try:
        registered = get_registered_tests(build_dir)
    except (
        OSError,
        ValueError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ) as exc:
        print(f"Cannot read CTest registry: {exc}", file=sys.stderr)
        return 2
    rows = []
    for case in cases:
        if args.preflight_only:
            status = (
                "REGISTERED_NOT_RUN" if case.name in registered else "NOT_REGISTERED"
            )
            rows.append(
                {\n                    **asdict(case),\n                    "status": status,\n                    "reason": "preflight-only: no execution",\n                }
            )
            print(f"{case.name}: {rows[-1]['status']}", flush=True)
            continue
        if case.name not in registered:
            rows.append(
                {
                    **asdict(case),
                    "status": "NOT_REGISTERED",
                    "reason": (
                        "CTest test absent; full validation requires"
                        " -DCFDX_ENABLE_LONG_VALIDATION=ON"
                    ),
                }
            )
        else:
            rows.append(run_case(case, build_dir, logs, args.jobs))
        print(f"{case.name}: {rows[-1]['status']}", flush=True)
    result = {
        "schema": "cfdx-qualification-matrix/v1",
        "utc": datetime.now(timezone.utc).isoformat(),
        "tier": args.tier,
        "build_dir": str(build_dir),
        "cases": rows,
        "all_pass": all(row["status"] == "PASS" for row in rows),
        "warning": "CTest pass is not global CFD physics qualification",
    }
    (logs / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return 0 if result["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
