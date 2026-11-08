"""JUnit evidence parser tests, independent of a C++ build."""
from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ctest_timings", ROOT / "scripts/ci/ctest_timings.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CTestTimingTests(unittest.TestCase):
    def test_summary_identifies_failure(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "junit.xml"
            p.write_text('<testsuite tests="2"><testcase name="one" time="0.123"/>'
                         '<testcase name="two" time="3.1"><failure/></testcase></testsuite>')
            result = MODULE.summarize([p])
            self.assertEqual(result["executed"], 2)
            self.assertEqual(result["failed_or_skipped"][0]["name"], "two")
            self.assertAlmostEqual(result["sum_test_seconds"], 3.223)

    def test_no_test_is_not_valid_evidence(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "empty.xml"
            p.write_text('<testsuite tests="0"/>')
            with self.assertRaisesRegex(ValueError, "no testcase"):
                MODULE.summarize([p])


if __name__ == "__main__":
    unittest.main()
