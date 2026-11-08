"""Offline regression tests for the fail-closed CTest lane planner."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "ctest_shard", ROOT / "scripts/ci/ctest_shard.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def fake(name, labels="long", **props):
    return {
        "name": name,
        "command": [name],
        "properties": [
            {"name": "LABELS", "value": labels.split(";")},
            *({"name": k, "value": v} for k, v in props.items()),
        ],
    }


class CTestShardTests(unittest.TestCase):
    def test_every_selected_name_once_deterministically(self):
        tests = [fake(f"test_{i}", TIMEOUT=i * 120) for i in range(1, 13)]
        tests.append(fake("quick_only", "unit"))
        a = MODULE.build_plan(tests, "long", 4, required_tests=("test_4",))
        b = MODULE.build_plan(
            list(reversed(tests)), "long", 4, required_tests=("test_4",)
        )
        self.assertEqual(a, b)
        assigned = [name for lane in a["lanes"] for name in lane["tests"]]
        self.assertCountEqual(assigned, [f"test_{i}" for i in range(1, 13)])
        self.assertNotIn("quick_only", assigned)

    def test_depends_is_not_split(self):
        tests = [
            fake("first", TIMEOUT=300),
            fake("second", DEPENDS=["first"]),
            fake("third", TIMEOUT=300),
            fake("fourth", TIMEOUT=300),
        ]
        plan = MODULE.build_plan(tests, "long", 2)
        self.assertTrue(
            any({"first", "second"}.issubset(set(l["tests"])) for l in plan["lanes"])
        )

    def test_fixture_setup_outside_label_is_allowed(self):
        tests = [
            fake("meshgen", "mesh", FIXTURES_SETUP=["generated"]),
            fake("a", DEPENDS=["meshgen"], FIXTURES_REQUIRED=["generated"]),
            fake("b", TIMEOUT=10),
        ]
        self.assertEqual(MODULE.build_plan(tests, "long", 2)["total_selected"], 2)

    def test_uncovered_dependency_is_rejected(self):
        tests = [fake("quick", "unit"), fake("a", DEPENDS=["quick"]), fake("b")]
        with self.assertRaisesRegex(ValueError, "outside long suite"):
            MODULE.build_plan(tests, "long", 2)

    def test_missing_required_test_is_failure(self):
        with self.assertRaisesRegex(ValueError, "required long test missing"):
            MODULE.build_plan([fake("a"), fake("b")], "long", 2, ("missing",))

    def test_disabled_gate_is_failure(self):
        with self.assertRaisesRegex(ValueError, "disabled long test"):
            MODULE.build_plan([fake("a", DISABLED=True), fake("b")], "long", 2)

    def test_zero_test_lane_is_failure(self):
        with self.assertRaisesRegex(ValueError, "only 1 tests"):
            MODULE.build_plan([fake("a")], "long", 4)

    def test_plan_writes_ctest_compatible_anchored_regex(self):
        plan = MODULE.build_plan([fake("test.a"), fake("test+b")], "long", 2)
        with tempfile.TemporaryDirectory() as td:
            MODULE.write_plan(plan, Path(td))
            for lane in plan["lanes"]:
                name = lane["tests"][0]
                regex = (
                    (Path(td) / f"long-shard-{lane['index']}.regex").read_text().strip()
                )
                import re
                self.assertTrue(re.fullmatch(regex, name))
                self.assertFalse(re.fullmatch(regex, "testaaaa"))

    def test_parse_measured_cost_file(self):
        with tempfile.TemporaryDirectory() as td:
            file = Path(td) / "CTestCostData.txt"
            file.write_text("a 4 3.5\nb 1 45\n---\na\n")
            self.assertEqual(MODULE.costs_from_ctest(file), {"a": 3.5, "b": 45})


if __name__ == "__main__":
    unittest.main()