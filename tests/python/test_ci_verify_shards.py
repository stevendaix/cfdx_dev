"""Regression tests for complete cross-runner long qualification coverage."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "verify_shards", ROOT / "scripts/ci/verify_shards.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def evidence(
    root: Path, failed: bool = False, duplicate: bool = False, omitted: bool = False
) -> Path:
    plan = {
        "suite_label": "long",
        "number_of_lanes": 2,
        "total_selected": 3,
        "lanes": [
            {"index": 0, "test_count": 2, "tests": ["a", "b"]},
            {"index": 1, "test_count": 1, "tests": ["c"]},
        ],
    }
    for shard, names in enumerate((["a", "b"], ["c"])):
        sub = root / f"cfdx-total-validation-long-evidence-{shard}" / "build-validation"
        (sub / "ci-shards").mkdir(parents=True)
        (sub / "ci-shards/long-manifest.json").write_text(json.dumps(plan))
        (sub / f"ci-shards/long-shard-{shard}.txt").write_text("\n".join(names) + "\n")
        cases = []
        for name in names:
            if omitted and name == "c":
                continue
            failure = "<failure/>" if failed and name == "b" else ""
            cases.append(f'<testcase name="{name}" time="1.0">{failure}</testcase>')
        if duplicate and shard == 1:
            cases.append('<testcase name="a" time="0.1"/>')
        (sub / f"ctest-long-{shard}.xml").write_text(
            "<testsuite>" + "".join(cases) + "</testsuite>"
        )
    return root


class VerifyShardsTests(unittest.TestCase):
    def test_complete_matrix_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = MODULE.verify(evidence(Path(tmp)), "long", 2)
            self.assertEqual(summary["long_tests_completed"], 3)

    def test_missing_shard_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = evidence(Path(tmp))
            import shutil

            shutil.rmtree(root / "cfdx-total-validation-long-evidence-1")
            with self.assertRaisesRegex(ValueError, "missing long-test artifact"):
                MODULE.verify(root, "long", 2)

    def test_skipped_test_fails(self):
        with (
            tempfile.TemporaryDirectory() as tmp,
            self.assertRaisesRegex(
                ValueError, "no CTest JUnit test cases|expected one execution"
            ),
        ):
            MODULE.verify(evidence(Path(tmp), omitted=True), "long", 2)

    def test_wrong_lane_is_rejected(self):
        with (
            tempfile.TemporaryDirectory() as tmp,
            self.assertRaisesRegex(ValueError, "wrong shard"),
        ):
            MODULE.verify(evidence(Path(tmp), duplicate=True), "long", 2)

    def test_failure_rejected(self):
        with (
            tempfile.TemporaryDirectory() as tmp,
            self.assertRaisesRegex(ValueError, "non-passing test"),
        ):
            MODULE.verify(evidence(Path(tmp), failed=True), "long", 2)


if __name__ == "__main__":
    unittest.main()
