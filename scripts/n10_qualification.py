#!/usr/bin/env python3
"""Run the reproducible N10 difficult-mesh qualification campaign.

This driver reuses the existing N10-labelled V&V executables. It does not
create a second smooth-field campaign and does not alter numerical tolerances.
The report records the exact CTest results and structured N10 diagnostics.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None
import time
from pathlib import Path

REQUIRED_TESTS = (
    "test_n10_difficult_mesh_robustness",
    "test_gradient_verification",
    "test_polyhedral_gradient_campaign",
    "test_nonorthogonal_skew_campaign",
    "test_nonorthogonal_laplacian_campaign",
    "test_polyhedral_laplacian_campaign",
)

_PREFIX = re.compile(r"^\s*Test\s+#\s*\d+\s*:\s*(\S+)")
_CTEST_PREFIX = re.compile(r"^\d+: ?")


def ctest(build_dir: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["ctest", *args],
        cwd=build_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def discover_tests(build_dir: Path) -> set[str]:
    result = ctest(build_dir, "-N")
    if result.returncode != 0:
        raise RuntimeError(f"ctest -N failed with exit code {result.returncode}\n{result.stdout}")
    return {
        match.group(1)
        for line in result.stdout.splitlines()
        if (match := _PREFIX.match(line))
    }


def run_test(build_dir: Path, name: str) -> dict[str, object]:
    start = time.monotonic()
    result = ctest(
        build_dir,
        "-V",
        "--output-on-failure",
        "--no-tests=error",
        "--tests-regex",
        f"^{name}$",
    )
    elapsed = time.monotonic() - start
    return {
        "name": name,
        "returncode": result.returncode,
        "status": "PASS" if result.returncode == 0 else "FAIL",
        "elapsed_s": round(elapsed, 3),
        "output": result.stdout,
    }


def parse_records(output: str, prefix: str) -> list[dict[str, object]]:
    pattern = re.compile(r"^" + re.escape(prefix) + r"\s+(.*)$")
    records: list[dict[str, object]] = []
    for raw in output.splitlines():
        line = _CTEST_PREFIX.sub("", raw, count=1)
        match = pattern.match(line)
        if not match:
            continue
        record: dict[str, object] = {}
        for token in match.group(1).split():
            key, sep, value = token.partition("=")
            if not sep:
                continue
            try:
                parsed: object = float(value) if any(c in value for c in ".eE") else int(value)
            except ValueError:
                parsed = value
            record[key] = parsed
        if record:
            records.append(record)
    return records


def extract_evidence(results: list[dict[str, object]]) -> dict[str, list[dict[str, object]]]:
    evidence: dict[str, list[dict[str, object]]] = {
        "quality": [],
        "solver": [],
        "polyhedral": [],
        "near_degenerate": [],
        "invalid": [],
    }
    prefixes = {
        "N10_QUALITY": "quality",
        "N10_SOLVER": "solver",
        "N10_POLYHEDRAL": "polyhedral",
        "N10_NEAR_DEGENERATE_VALID": "near_degenerate",
        "N10_INVALID": "invalid",
    }
    for result in results:
        output = str(result["output"])
        for prefix, category in prefixes.items():
            evidence[category].extend(parse_records(output, prefix))
    return evidence


def audit_evidence(evidence: dict[str, list[dict[str, object]]]) -> dict[str, object]:
    required = {
        "quality": ("case", "cells", "max_skewness", "max_nonorth_deg", "max_aspect", "min_volume", "max_volume"),
        "solver": ("case", "cells", "iterations", "reported_relative_residual", "true_relative_residual"),
        "polyhedral": ("cells", "max_skewness", "max_nonorth_deg", "max_aspect", "linear_gradient_Linf"),
        "near_degenerate": ("max_aspect", "min_volume"),
        "invalid": ("rejected_errors",),
    }
    missing: dict[str, list[int]] = {}
    checked = 0
    complete = 0
    for category, fields in required.items():
        records = evidence.get(category, [])
        for index, record in enumerate(records):
            checked += 1
            absent = [field for field in fields if field not in record]
            if absent:
                missing[f"{category}[{index}]"] = absent
            else:
                complete += 1
    empty = [category for category, records in evidence.items() if not records]
    return {
        "status": "INCOMPLETE" if missing or empty else "COMPLETE",
        "records_checked": checked,
        "records_complete": complete,
        "missing_fields": missing,
        "empty_categories": empty,
        "policy": "diagnostic_only",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--fixture-manifest", type=Path, default=None)
    args = parser.parse_args()

    build_dir = args.build_dir.resolve()
    if not (build_dir / "CTestTestfile.cmake").exists():
        print(f"error: {build_dir} is not a configured CMake build directory", file=sys.stderr)
        return 2

    report_path = args.report.resolve() if args.report else build_dir / "n10_qualification.json"
    repo_root = Path(__file__).resolve().parents[1]
    manifest_path = (
        args.fixture_manifest
        if args.fixture_manifest is not None
        else repo_root / "tests/fixtures/mesh_sources.yaml"
    ).resolve()
    if not manifest_path.exists():
        print(f"error: fixture manifest does not exist: {manifest_path}", file=sys.stderr)
        return 2
    if yaml is None:
        print("error: PyYAML is required to parse the N10 fixture manifest", file=sys.stderr)
        return 2
    try:
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"error: unable to parse fixture manifest: {exc}", file=sys.stderr)
        return 2

    if not isinstance(manifest, dict):
        print("error: fixture manifest root must be a mapping", file=sys.stderr)
        return 2
    fixtures = manifest.get("fixtures")
    if not isinstance(fixtures, list):
        print("error: fixture manifest 'fixtures' must be a list", file=sys.stderr)
        return 2

    verified_fixtures: list[dict[str, object]] = []
    for index, item in enumerate(fixtures):
        if not isinstance(item, dict):
            print(f"error: fixture manifest entry {index} must be a mapping", file=sys.stderr)
            return 2
        if item.get("status") != "verified_reference":
            continue
        required_fields = ("id", "source", "acquisition", "sha256")
        missing_fields = [field for field in required_fields if not item.get(field)]
        if missing_fields:
            print(
                f"error: verified fixture entry {index} is missing required fields: "
                + ", ".join(missing_fields),
                file=sys.stderr,
            )
            return 2
        verified_fixtures.append(
            {
                "id": item["id"],
                "source": item["source"],
                "acquisition": item["acquisition"],
                "sha256": item["sha256"],
                "status": item["status"],
            }
        )
    if not verified_fixtures:
        print("error: fixture manifest contains no verified_reference fixtures", file=sys.stderr)
        return 2

    available = discover_tests(build_dir)
    missing = [name for name in REQUIRED_TESTS if name not in available]
    if missing:
        report = {
            "campaign": "N10 difficult-mesh qualification",
            "status": "INCOMPLETE",
            "required_tests": list(REQUIRED_TESTS),
            "missing_tests": missing,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 2

    results: list[dict[str, object]] = []
    for name in REQUIRED_TESTS:
        print(f"\n=== N10 qualification: {name} ===", flush=True)
        result = run_test(build_dir, name)
        results.append(result)
        print(result["output"], end="", flush=True)
        if result["returncode"] != 0:
            break

    evidence = extract_evidence(results)
    coverage = audit_evidence(evidence)
    failed = [r["name"] for r in results if r["status"] == "FAIL"]
    complete = len(results) == len(REQUIRED_TESTS)
    status = "PASS" if complete and not failed and coverage["status"] == "COMPLETE" else "FAIL"

    try:
        git_revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        git_revision = "UNKNOWN"

    report = {
        "campaign": "N10 difficult-mesh qualification",
        "git_revision": git_revision,
        "status": status,
        "required_tests": list(REQUIRED_TESTS),
        "completed_tests": len(results),
        "failed_tests": failed,
        "results": results,
        "evidence": evidence,
        "evidence_coverage": coverage,
        "scope": {
            "geometry_robustness": True,
            "reconstruction_and_pde_accuracy": "reused from N2/N3 campaigns",
            "solver_robustness": "reported by dedicated N10 campaign",
            "imported_production_fixtures": "verified_reference_manifest_integrated",
            "verified_fixture_count": len(verified_fixtures),
            "verified_fixtures": verified_fixtures,
        },
        "policy": {
            "reuses_existing_vv_tests": True,
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_fallbacks": False,
            "universal_polyhedral_accuracy_claim": False,
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nN10 qualification report: {report_path}")
    print(f"N10 qualification status: {status}")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
