"""Do not turn partial qualification tracking into a false model PASS."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/qualification_portfolio.py"
SPEC = importlib.util.spec_from_file_location("qualification_portfolio", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
import sys

sys.path.append(str(ROOT / "scripts"))
portfolio = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(portfolio)

REGISTRY = portfolio.read_json(portfolio.REGISTRY)
LINKS = portfolio.read_json(portfolio.LINKS)


class QualificationPortfolioTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.manifest = {
            n for info in LINKS["cases"].values()
            for n in info["contract_tests"] + info["campaign_tests"]
        } | {
            n for info in LINKS["unregistered_backlog"]
            for n in info["existing_checks"]
        }

    def test_exact_25_case_mapping_and_backlog_are_distinct(self):
        portfolio.validate_links(REGISTRY, LINKS)
        self.assertEqual(len(LINKS["cases"]), 25)
        self.assertGreaterEqual(len(LINKS["unregistered_backlog"]), 5)
        self.assertEqual(set(LINKS["cases"]), {c["id"] for c in REGISTRY["cases"]})

    def test_duplicate_or_missing_case_mapping_is_rejected(self):
        alternate = json.loads(json.dumps(LINKS))
        alternate["cases"].pop("INT-CHANNEL")
        with self.assertRaisesRegex(portfolio.PortfolioError, "case mapping drift"):
            portfolio.validate_links(REGISTRY, alternate)
        alternate = json.loads(json.dumps(LINKS))
        alternate["cases"]["LAM-COUETTE"]["campaign_tests"].append("test_couette_quick")
        with self.assertRaisesRegex(portfolio.PortfolioError, "duplicate contract/campaign"):
            portfolio.validate_links(REGISTRY, alternate)

    def test_manifest_and_junit_duplicate_cannot_become_green(self):
        manifest = self.base / "manifest.json"
        manifest.write_text(json.dumps({"tests": [{"name": n} for n in sorted(self.manifest)]}))
        self.assertEqual(portfolio.read_manifest(manifest), {name: True for name in self.manifest})
        junit = self.base / "ctest.xml"
        junit.write_text('<testsuite tests="2"><testcase name="test_couette_quick"/><testcase name="test_couette_quick"/></testsuite>')
        self.assertEqual(portfolio.read_junits([junit])["test_couette_quick"], "DUPLICATE")
        self.assertEqual(portfolio.aggregate_status({"test_couette_quick": "DUPLICATE"}), "AMBIGUOUS")

    def test_disabled_ctest_manifest_anchor_is_not_a_pass(self):
        manifest = self.base / "manifest.json"
        manifest.write_text(json.dumps({"tests": [{"name": "test_couette_quick", "properties": [{"name": "DISABLED", "value": True}]}]}))
        configured = portfolio.read_manifest(manifest)
        self.assertEqual(portfolio.test_state("test_couette_quick", configured, {}), "DISABLED")
        self.assertEqual(portfolio.aggregate_status({"test_couette_quick": "DISABLED"}), "DISABLED")

    def test_skipped_fail_unregistered_and_not_run_are_not_pass(self):
        junit = self.base / "ctest.xml"
        junit.write_text('<testsuite><testcase name="test_couette_quick"><skipped/></testcase><testcase name="test_poiseuille_quick"><failure/></testcase></testsuite>')
        statuses = portfolio.read_junits([junit])
        self.assertEqual(statuses["test_couette_quick"], "SKIPPED")
        self.assertEqual(statuses["test_poiseuille_quick"], "FAIL")
        self.assertEqual(portfolio.test_state("missing", self.manifest, statuses), "NOT_REGISTERED")
        self.assertEqual(portfolio.test_state("test_channel_quick", self.manifest, statuses), "NOT_RUN")
        self.assertNotEqual(portfolio.aggregate_status({"x": "SKIPPED"}), "PASS")
        self.assertNotEqual(portfolio.aggregate_status({"x": "FAIL"}), "PASS")

    def test_historical_pass_and_green_tests_never_promote_physics(self):
        with_all_pass = {name: "PASS" for name in self.manifest}
        board, evidence = portfolio.make_portfolio(
            REGISTRY, LINKS, self.manifest, with_all_pass, self.base / "nonexistent_evidence"
        )
        self.assertEqual(board["summary"]["physical_evidence_gates_pass"], 0)
        self.assertEqual(evidence["coverage"]["not_qualified"], 25)
        for case in board["cases"]:
            self.assertEqual(case["physical_evidence_status"], "NOT_QUALIFIED")
        ghia = next(c for c in board["cases"] if c["id"] == "INC-GHIA")
        self.assertEqual(ghia["registry_status"], "PASS")
        self.assertEqual(ghia["campaign_status"], "PASS")
        self.assertNotEqual(ghia["physical_evidence_status"], "EVIDENCE_GATES_PASS")

    def test_no_junit_means_tests_not_run(self):
        board, _ = portfolio.make_portfolio(
            REGISTRY, LINKS, self.manifest, {}, self.base / "missing_evidence"
        )
        self.assertEqual(board["summary"]["provided_junit_test_results"], 0)
        self.assertEqual(board["summary"]["registered_missing_ctest_anchors"], [])
        self.assertEqual(next(c for c in board["cases"] if c["id"] == "INC-GHIA")["campaign_status"], "NOT_RUN")
        self.assertEqual(next(c for c in board["cases"] if c["id"] == "RAD-S2S")["contract_status"], "NOT_RUN")

    def test_markdown_diff_records_changes_and_next_actions(self):
        before, _ = portfolio.make_portfolio(REGISTRY, LINKS, self.manifest, {}, self.base / "none")
        after, _ = portfolio.make_portfolio(REGISTRY, LINKS, self.manifest, {"test_couette_quick": "PASS"}, self.base / "none")
        md = portfolio.markdown_report(after, before)
        self.assertIn("Changements depuis le précédent bilan", md)
        self.assertIn("`LAM-COUETTE`", md)
        self.assertIn("NOT_RUN", md)
        self.assertIn("PASS", md)
        self.assertIn("TH-RAYLEIGH-BENARD", md)
        self.assertIn("0/25", md)

    def test_cli_writes_three_machine_readable_outputs(self):
        manifest = self.base / "manifest.json"
        manifest.write_text(json.dumps({"tests": [{"name": n} for n in sorted(self.manifest)]}))
        output = self.base / "dashboard"
        code = portfolio.main(["--ctest-manifest", str(manifest), "--evidence-dir", str(self.base / "none"), "--output-dir", str(output), "--strict-links"])
        self.assertEqual(code, 0)
        for filename in ("qualification_portfolio.md", "qualification_portfolio.csv", "qualification_portfolio.json", "physical_evidence_gate.json"):
            self.assertTrue((output / filename).is_file())
        report = json.loads((output / "qualification_portfolio.json").read_text())
        self.assertEqual(len(report["cases"]), 25)
        code = portfolio.main(["--ctest-manifest", str(manifest), "--evidence-dir", str(self.base / "none"), "--output-dir", str(output), "--require-evidence"])
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
