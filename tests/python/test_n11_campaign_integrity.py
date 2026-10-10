#!/usr/bin/env python3
"""Regression tests for N11 exact-HEAD campaign evidence integrity."""

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from n11_conservation_boundedness import expected_sha_matches


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
