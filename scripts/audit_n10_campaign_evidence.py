"""Fail closed unless the complete N10 campaign report matches the tested HEAD."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

REQUIRED_TESTS = {
    "test_n10_difficult_mesh_robustness",
    "test_gradient_verification",
    "test_polyhedral_gradient_campaign",
    "test_nonorthogonal_skew_campaign",
    "test_nonorthogonal_laplacian_campaign",
    "test_polyhedral_laplacian_campaign",
}


def audit(report: object, expected_revision: str) -> list[str]:
    errors: list[str] = []
    if not isinstance(report, dict):
        return ["campaign report root must be a JSON object"]

    revision = report.get("git_revision")
    if not isinstance(revision, str) or not re.fullmatch(
        r"[0-9a-fA-F]{40,64}", revision
    ):
        errors.append(f"git_revision is missing or malformed: {revision!r}")
    elif revision.lower() != expected_revision.lower():
        errors.append(
            f"git_revision={revision} does not match tested HEAD={expected_revision}"
        )

    if report.get("status") != "PASS":
        errors.append(f"campaign status={report.get('status')!r}, expected 'PASS'")

    required_tests = report.get("required_tests")
    if (
        not isinstance(required_tests, list)
        or any(not isinstance(name, str) for name in required_tests)
        or set(required_tests) != REQUIRED_TESTS
    ):
        errors.append("required_tests does not match the authoritative N10 test set")
    if isinstance(required_tests, list) and all(
        isinstance(name, str) for name in required_tests
    ) and len(required_tests) != len(set(required_tests)):
        errors.append("required_tests contains duplicate entries")

    completed = report.get("completed_tests")
    if completed != len(REQUIRED_TESTS):
        errors.append(f"completed_tests={completed!r}, expected {len(REQUIRED_TESTS)}")

    results = report.get("results")
    if not isinstance(results, list):
        errors.append("results is missing or is not a list")
    else:
        result_names = [
            row.get("name") for row in results if isinstance(row, dict)
        ]
        if len(result_names) != len(results):
            errors.append("results contains a non-object entry")
        if (
            any(not isinstance(name, str) for name in result_names)
            or set(result_names) != REQUIRED_TESTS
            or len(result_names) != len(REQUIRED_TESTS)
        ):
            errors.append("results do not contain each required N10 test exactly once")
        for row in results:
            if isinstance(row, dict) and (
                row.get("status") != "PASS" or row.get("returncode") != 0
            ):
                errors.append(
                    f"test {row.get('name')!r} status={row.get('status')!r}, "
                    f"returncode={row.get('returncode')!r}"
                )

    if report.get("failed_tests") != []:
        errors.append(f"failed_tests={report.get('failed_tests')!r}, expected []")

    coverage = report.get("evidence_coverage")
    if not isinstance(coverage, dict) or coverage.get("status") != "COMPLETE":
        errors.append("evidence_coverage is missing or not COMPLETE")

    policy = report.get("policy")
    if not isinstance(policy, dict):
        errors.append("campaign policy is missing")
    else:
        for key in (
            "changes_numerical_tolerances",
            "disables_validation",
            "silent_fallbacks",
        ):
            if policy.get(key) is not False:
                errors.append(f"policy {key}={policy.get(key)!r}, expected false")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--expected-revision", required=True)
    args = parser.parse_args()

    if not re.fullmatch(r"[0-9a-fA-F]{40,64}", args.expected_revision):
        print(
            "error: --expected-revision must be a full Git commit SHA",
            file=sys.stderr,
        )
        return 2
    try:
        report = json.loads(args.report.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(
            f"error: cannot read complete N10 campaign report: {exc}",
            file=sys.stderr,
        )
        return 2

    errors = audit(report, args.expected_revision)
    if errors:
        print("N10 exact-HEAD campaign evidence: FAIL")
        for error in errors:
            print(f"- {error}")
        return 1

    print(f"N10 exact-HEAD campaign evidence: PASS ({args.expected_revision})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
