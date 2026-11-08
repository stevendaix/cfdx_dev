"""Workflow integrity checks with no dependency on the CFDX C++ solver."""
from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ci_policy", ROOT / "scripts/ci/ci_policy.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CIPolicyTests(unittest.TestCase):
    def test_all_workflow_definitions_parse_and_partition_gate_complete(self):
        result = MODULE.audit(ROOT / ".github/workflows")
        self.assertGreaterEqual(result["workflows_checked"], 20)

    def test_duplicate_yaml_key_is_a_fatal_error(self):
        with tempfile.TemporaryDirectory() as td:
            file = Path(td) / "bad.yml"
            file.write_text('name: hi\non:\n  push:\n  push:\njobs:\n  foo: {}\n')
            with self.assertRaisesRegex(ValueError, "duplicate YAML key"):
                MODULE.load(file)

    def test_github_on_key_is_not_automatic_yaml_boolean(self):
        with tempfile.TemporaryDirectory() as td:
            file = Path(td) / "events.yml"
            file.write_text('name: hi\non:\n  workflow_dispatch:\njobs:\n  foo: {}\n')
            workflow = MODULE.load(file)
            self.assertIn("on", workflow)


if __name__ == "__main__":
    unittest.main()