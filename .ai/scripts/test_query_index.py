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
            """)
            db.execute("INSERT INTO index_metadata VALUES ('repository_revision','unknown')")
            db.execute("INSERT INTO files VALUES (1,'a.py',?)", (MODULE.sha256(source),))
            db.execute("INSERT INTO symbols VALUES (1,1,'function','alpha','alpha',1,4)")
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


if __name__ == "__main__":
    unittest.main()
