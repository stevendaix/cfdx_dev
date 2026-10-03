#!/usr/bin/env python3
"""Tests for deterministic code-index validation and queries."""
from __future__ import annotations

import importlib.util
import pathlib
import sqlite3
import tempfile
import unittest


SCRIPT = pathlib.Path(__file__).with_name("query_index.py")
SPEC = importlib.util.spec_from_file_location("query_index", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class QueryIndexTests(unittest.TestCase):
    def _db(self, root: pathlib.Path) -> pathlib.Path:
        source = root / "a.py"
        source.write_text("def alpha():\n    return 1\n", encoding="utf-8")
        db_path = root / "index.sqlite"
        with sqlite3.connect(db_path) as db:
            db.executescript("""
                CREATE TABLE index_metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
                CREATE TABLE files(id INTEGER PRIMARY KEY, path TEXT, content_hash TEXT);
                CREATE TABLE symbols(id INTEGER PRIMARY KEY, file_id INTEGER, kind TEXT,
                    qualified_name TEXT, name TEXT, line_start INTEGER, column_start INTEGER);
                CREATE TABLE code_references(id INTEGER PRIMARY KEY, source_symbol_id INTEGER,
                    target_symbol_id INTEGER, file_id INTEGER, kind TEXT, line INTEGER, column INTEGER);
                CREATE TABLE dependencies(id INTEGER PRIMARY KEY, source_file_id INTEGER,
                    target_path TEXT, kind TEXT);
                CREATE TABLE tests(id INTEGER PRIMARY KEY, file_id INTEGER, name TEXT,
                    framework TEXT, line_start INTEGER, line_end INTEGER);
                CREATE TABLE validation_cases(id INTEGER PRIMARY KEY, path TEXT, name TEXT,
                    category TEXT, status TEXT, source_revision TEXT);
                CREATE TABLE test_links(test_id INTEGER, symbol_id INTEGER, relation TEXT);
                CREATE TABLE validation_links(validation_case_id INTEGER, symbol_id INTEGER, relation TEXT);
            """)
            db.execute("INSERT INTO index_metadata VALUES ('repository_revision','unknown')")
            db.execute("INSERT INTO files VALUES (1,'a.py',?)", (MODULE.sha256(source),))
            db.execute("INSERT INTO symbols VALUES (1,1,'function','alpha','alpha',1,4)")
            db.execute("INSERT INTO tests VALUES (1,1,'alpha_test','python',1,3)")
            db.execute("INSERT INTO validation_cases VALUES (1,'tests/validation/alpha','Alpha case','unit','verified','unknown')")
            db.execute("INSERT INTO test_links VALUES (1,1,'defines-symbol')")
            db.execute("INSERT INTO validation_links VALUES (1,1,'explicit-test-reference')")
            db.commit()
        return db_path

    def test_validate_accepts_unchanged_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            db = self._db(root)
            with sqlite3.connect(db) as conn:
                self.assertEqual(MODULE.validate(conn, root), [])

    def test_validate_detects_changed_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            db = self._db(root)
            (root / "a.py").write_text("def changed():\n    return 2\n", encoding="utf-8")
            with sqlite3.connect(db) as conn:
                errors = MODULE.validate(conn, root)
            self.assertIn("indexed file changed: a.py", errors)

    def test_symbol_query(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            db = self._db(root)
            with sqlite3.connect(db) as conn:
                rows = MODULE.symbols(conn, "alp")
            self.assertEqual(rows[0][1], "alpha")

    def test_test_query_reports_explicit_symbol_link(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            db = self._db(root)
            with sqlite3.connect(db) as conn:
                rows = MODULE.tests(conn, "alpha")
            self.assertEqual(rows[0][0:3], ("alpha_test", "python", "a.py"))
            self.assertEqual(rows[0][4], "alpha")
            self.assertEqual(rows[0][5], "defines-symbol")

    def test_validation_query_reports_explicit_symbol_link(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            db = self._db(root)
            with sqlite3.connect(db) as conn:
                rows = MODULE.validations(conn, "Alpha")
            self.assertEqual(rows[0][0:4], ("Alpha case", "tests/validation/alpha", "unit", "verified"))
            self.assertEqual(rows[0][5], "alpha")
            self.assertEqual(rows[0][6], "explicit-test-reference")


if __name__ == "__main__":
    unittest.main()
