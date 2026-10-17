"""Regression tests for N11 exact-HEAD campaign evidence integrity."""

import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "n11_conservation_boundedness.py"
SPEC = importlib.util.spec_from_file_location("n11_conservation_boundedness", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
CAMPAIGN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CAMPAIGN)
expected_sha_matches = CAMPAIGN.expected_sha_matches


class N11CampaignIntegrityTests(unittest.TestCase):
    def test_matching_full_sha_is_accepted(self):
        sha = "a" * 40
        self.assertTrue(expected_sha_matches(sha, sha))

    def test_different_sha_is_rejected(self):
        self.assertFalse(expected_sha_matches("a" * 40, "b" * 40))

    def test_missing_expected_sha_allows_local_run(self):
        self.assertTrue(expected_sha_matches("a" * 40, None))

    def test_empty_expected_sha_allows_local_run(self):
        self.assertTrue(expected_sha_matches("a" * 40, ""))

class ParseQuantitativeRecordsTests(unittest.TestCase):
    def test_plain_record_line_is_parsed(self):
        records = CAMPAIGN.parse_quantitative_records(
            "MODEL_RESULT case=simple_algorithm observed_order=2.0"
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["record_type"], "MODEL_RESULT")
        self.assertEqual(records[0]["case"], "simple_algorithm")
        self.assertEqual(records[0]["observed_order"], 2.0)

    def test_ctest_verbose_prefix_is_stripped(self):
        records = CAMPAIGN.parse_quantitative_records(
            "106: POISEUILLE_RESULT level=8 l2=1.25e-3"
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["record_type"], "POISEUILLE_RESULT")
        self.assertEqual(records[0]["level"], 8.0)
        self.assertEqual(records[0]["l2"], 1.25e-3)

    def test_prefixed_non_record_lines_are_ignored(self):
        records = CAMPAIGN.parse_quantitative_records(
            "24: PASS: test_case\n25: ITER 1 continuity=1e-12"
        )
        self.assertEqual(records, [])

    def test_prefix_stripping_preserves_record_content(self):
        plain = CAMPAIGN.parse_quantitative_records("GHIA re=100 max=0.05")
        prefixed = CAMPAIGN.parse_quantitative_records("7: GHIA re=100 max=0.05")
        self.assertEqual(plain, prefixed)

if __name__ == "__main__":
    unittest.main()

