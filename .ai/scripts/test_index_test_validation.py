#!/usr/bin/env python3
"""Tests for conservative test/validation inventory discovery."""
from __future__ import annotations

import importlib.util
import pathlib
import sqlite3
import sys
import tempfile
import unittest

SCRIPT = pathlib.Path(__file__).with_name("index_test_validation.py")
SPEC = importlib.util.spec_from_file_location("index_test_validation", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class TestValidationIndexerTests(unittest.TestCase):
    def test_discovery_is_conservative_and_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "tests/unit").mkdir(parents=True)
            (root / "tests/validation/couette").mkdir(parents=True)
            (root / "tests/unit/test_solver.cpp").write_text("", encoding="utf-8")
            (root / "tests/unit/solver.cpp").write_text("", encoding="utf-8")
            (root / "tests/validation/couette/case.yaml").write_text("", encoding="utf-8")

            self.assertEqual(
                MODULE.discover_tests(root),
                [("tests/unit/test_solver.cpp", "source-discovery")],
            )
            self.assertEqual(
                MODULE.discover_validation_cases(root),
                [("tests/validation/couette", "couette", "directory-discovery")],
            )

    def test_populates_existing_schema_without_inventing_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "tests").mkdir()
            (root / "tests/validation/couette").mkdir(parents=True)
            test_file = root / "tests/test_solver.py"
            test_file.write_text("", encoding="utf-8")

            db_path = root / "index.sqlite"
            schema = pathlib.Path(__file__).parents[1] / "index" / "schema.sql"
            with sqlite3.connect(db_path) as db:
                db.executescript(schema.read_text(encoding="utf-8"))
                db.execute(
                    "INSERT INTO files(path, language, content_hash, size_bytes, revision, indexed_at) "
                    "VALUES (?, ?, ?, ?, ?, ?)",
                    ("tests/test_solver.py", "python", "hash", 0, "rev", "now"),
                )

            counts = MODULE.index_test_validation(root, db_path)
            self.assertEqual(counts, (1, 1))
            with sqlite3.connect(db_path) as db:
                test = db.execute("SELECT name, framework FROM tests").fetchone()
                case = db.execute("SELECT path, name, status FROM validation_cases").fetchone()
                links = db.execute("SELECT COUNT(*) FROM test_links").fetchone()[0]
            self.assertEqual(test, ("test_solver", "source-discovery"))
            self.assertEqual(case, ("tests/validation/couette", "couette", "discovered"))
            self.assertEqual(links, 0)


if __name__ == "__main__":
    unittest.main()
