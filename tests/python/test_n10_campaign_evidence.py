"""Regression tests for the fail-closed N10 campaign evidence gate."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/audit_n10_campaign_evidence.py"
SPEC = importlib.util.spec_from_file_location("audit_n10_campaign_evidence", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

SHA = "a" * 40
TESTS = sorted(MODULE.REQUIRED_TESTS)


def valid_report() -> dict[str, object]:
    return {
        "git_revision": SHA,
        "status": "PASS",
        "required_tests": TESTS.copy(),
        "completed_tests": len(TESTS),
        "failed_tests": [],
        "results": [{"name": name, "status": "PASS", "returncode": 0} for name in TESTS],
        "evidence_coverage": {"status": "COMPLETE"},
        "policy": {
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_fallbacks": False,
        },
    }


class CampaignEvidenceTests(unittest.TestCase):
    def test_accepts_complete_report_for_exact_head(self) -> None:
        self.assertEqual(MODULE.audit(valid_report(), SHA), [])

    def test_rejects_wrong_head(self) -> None:
        errors = MODULE.audit(valid_report(), "b" * 40)
        self.assertTrue(any("does not match tested HEAD" in error for error in errors))

    def test_rejects_failed_test_even_if_report_claims_pass(self) -> None:
        report = valid_report()
        report["results"][0]["status"] = "FAIL"  # type: ignore[index]
        errors = MODULE.audit(report, SHA)
        self.assertTrue(any("status=" in error for error in errors))

    def test_rejects_incomplete_evidence(self) -> None:
        report = valid_report()
        report["evidence_coverage"] = {"status": "INCOMPLETE"}
        errors = MODULE.audit(report, SHA)
        self.assertTrue(any("not COMPLETE" in error for error in errors))

    def test_rejects_missing_tests_and_policy_regressions(self) -> None:
        report = valid_report()
        report["required_tests"] = TESTS[:-1]
        report["policy"]["silent_fallbacks"] = True  # type: ignore[index]
        errors = MODULE.audit(report, SHA)
        self.assertTrue(any("authoritative N10 test set" in error for error in errors))
        self.assertTrue(any("silent_fallbacks" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
