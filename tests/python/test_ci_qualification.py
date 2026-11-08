"""Fail-closed tests for the activation path of targeted CFDX qualification."""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "qualification", ROOT / "scripts/ci/qualification.py"
)
assert SPEC and SPEC.loader
Q = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(Q)
CATALOG = Q.read_catalog()


def make_test_case(name: str, **properties):
    return {
        "name": name,
        "properties": [{"name": k, "value": v} for k, v in properties.items()],
    }


class ProfileTests(unittest.TestCase):
    def test_catalog_has_required_domains_and_unique_labels(self):
        self.assertEqual(
            set(CATALOG),
            {
                "smoke",
                "naca",
                "momentum",
                "sst",
                "thermal",
                "wall-distance",
                "linalg",
                "n9",
                "ghia",
            },
        )
        self.assertEqual(len({p["label"] for p in CATALOG.values()}), len(CATALOG))

    def test_quick_profile_is_exact_not_fuzzy(self):
        fake = [make_test_case(n) for n in CATALOG["smoke"]["quick"]]
        fake.append(make_test_case("test_huge_campaign", LABELS=["long"]))
        plan = Q.plan_profile(fake, "smoke", "quick", CATALOG)
        self.assertEqual(plan["test_count"], len(CATALOG["smoke"]["quick"]))
        self.assertNotIn("test_huge_campaign", plan["selected_tests"])

    def test_mandatory_anchor_missing_rejects_plan(self):
        selected = [
            make_test_case(n)
            for n in CATALOG["naca"]["quick"]
            if n != "test_naca0012_quick"
        ]
        with self.assertRaisesRegex(
            Q.QualificationError, "required qualification tests missing"
        ):
            Q.plan_profile(selected, "naca", "quick", CATALOG)

    def test_extended_requires_hard_physical_anchor(self):
        selected = [make_test_case(n) for n in CATALOG["n9"]["quick"]]
        with self.assertRaisesRegex(Q.QualificationError, "does not select mandatory"):
            Q.plan_profile(selected, "n9", "extended", CATALOG)

    def test_mandatory_long_not_selected_by_quick(self):
        fake = [make_test_case(n) for n in CATALOG["ghia"]["quick"]]
        fake[0]["properties"].append({"name": "LABELS", "value": ["long"]})
        with self.assertRaisesRegex(Q.QualificationError, "quick profile.*long test"):
            Q.plan_profile(fake, "ghia", "quick", CATALOG)

    def test_disabled_selected_test_is_not_success(self):
        fake = [make_test_case(n) for n in CATALOG["smoke"]["quick"]]
        fake[1]["properties"].append({"name": "DISABLED", "value": True})
        with self.assertRaisesRegex(Q.QualificationError, "disabled qualification"):
            Q.plan_profile(fake, "smoke", "quick", CATALOG)

    def test_ctest_dependencies_and_fixtures_are_included_once(self):
        fake = [make_test_case(n) for n in CATALOG["smoke"]["quick"]]
        fake[0]["properties"].extend(
            [
                {"name": "DEPENDS", "value": ["prereq"]},
                {"name": "FIXTURES_REQUIRED", "value": ["mesh"]},
            ]
        )
        fake += [
            make_test_case("prereq"),
            make_test_case("mesh_setup", FIXTURES_SETUP=["mesh"]),
        ]
        plan = Q.plan_profile(fake, "smoke", "quick", CATALOG)
        self.assertEqual(plan["prerequisites"], ["mesh_setup", "prereq"])
        self.assertEqual(len(plan["selected_tests"]), len(set(plan["selected_tests"])))

    def test_missing_fixture_is_fatal(self):
        fake = [make_test_case(n) for n in CATALOG["smoke"]["quick"]]
        fake[0]["properties"].append(
            {"name": "FIXTURES_REQUIRED", "value": ["missing"]}
        )
        with self.assertRaisesRegex(Q.QualificationError, "missing fixture producer"):
            Q.plan_profile(fake, "smoke", "quick", CATALOG)

    def test_escaped_regex_and_saved_plan(self):
        plan = {"selected_tests": ["test_special_[1]", "ordinary"], "profile": "smoke"}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            Q.save_plan(plan, path)
            pattern = re.compile((path / "tests.regex").read_text().strip())
            self.assertTrue(pattern.match("test_special_[1]"))
            self.assertFalse(pattern.match("test_special_1"))
            self.assertTrue(pattern.match("ordinary"))
            self.assertFalse(pattern.match("unrelated"))

    def test_fully_empty_manifest_is_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            path.write_text(json.dumps({"tests": []}))
            with self.assertRaisesRegex(Q.QualificationError, "absent/empty"):
                Q.read_manifest(Path(tmp), path)

    def test_bad_profile_scope_rejected(self):
        fake = [make_test_case(n) for n in CATALOG["smoke"]["quick"]]
        with self.assertRaisesRegex(Q.QualificationError, "unknown qualification"):
            Q.plan_profile(fake, "does-not-exist", "quick", CATALOG)
        with self.assertRaisesRegex(
            Q.QualificationError, "unknown qualification scope"
        ):
            Q.plan_profile(fake, "smoke", "impossible", CATALOG)

    def test_missing_junit_is_fatal_for_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path / "tests.txt").write_text("a\nb\n")
            junit = path / "results.xml"
            junit.write_text('<testsuite><testcase name="a" time="1"/></testsuite>')
            cmd = [
                sys.executable,
                str(ROOT / "scripts/ci/ctest_timings.py"),
                str(junit),
                "--require-names-file",
                str(path / "tests.txt"),
            ]
            self.assertNotEqual(
                subprocess.run(cmd, check=False, capture_output=True).returncode, 0
            )
            junit.write_text(
                '<testsuite><testcase name="a" time="1"/><testcase name="b" time="2"/></testsuite>'
            )
            self.assertEqual(
                subprocess.run(cmd, check=False, capture_output=True).returncode, 0
            )
            junit.write_text(
                '<testsuite><testcase name="a" time="1"><skipped/></testcase><testcase name="b" time="2"/></testsuite>'
            )
            self.assertNotEqual(
                subprocess.run(cmd, check=False, capture_output=True).returncode, 0
            )


class TriggerTests(unittest.TestCase):
    def pr(self, labels):
        return {"pull_request": {"labels": [{"name": l} for l in labels]}}

    def test_manual_quick_and_total(self):
        self.assertEqual(
            Q.resolve("workflow_dispatch", {}, "naca", "quick", CATALOG)[:2],
            ("naca", "quick"),
        )
        self.assertEqual(
            Q.resolve("workflow_dispatch", {}, "total", "extended", CATALOG)[:2],
            ("total", "extended"),
        )

    def test_pr_label_runs_extended_ctest(self):
        self.assertEqual(
            Q.resolve("pull_request", self.pr(["qualification:naca"]), "", "", CATALOG)[
                :2
            ],
            ("naca", "extended"),
        )

    def test_no_pr_label_is_not_fake_green_qualification(self):
        self.assertEqual(
            Q.resolve("pull_request", self.pr(["bug"]), "", "", CATALOG)[0], "none"
        )

    def test_total_label_is_not_double_executed(self):
        self.assertEqual(
            Q.resolve(
                "pull_request",
                self.pr(["qualification:naca", "validation-total"]),
                "",
                "",
                CATALOG,
            )[0],
            "none",
        )

    def test_multiple_pr_labels_fail_closed(self):
        with self.assertRaisesRegex(
            Q.QualificationError, "multiple qualification labels"
        ):
            Q.resolve(
                "pull_request",
                self.pr(["qualification:naca", "qualification:thermal"]),
                "",
                "",
                CATALOG,
            )

    def test_unknown_manual_selection_fails_closed(self):
        with self.assertRaisesRegex(Q.QualificationError, "unknown manual profile"):
            Q.resolve("workflow_dispatch", {}, "not-real", "quick", CATALOG)


if __name__ == "__main__":
    unittest.main()
