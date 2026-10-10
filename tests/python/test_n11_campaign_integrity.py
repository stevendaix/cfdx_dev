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


if __name__ == "__main__":
    unittest.main()
