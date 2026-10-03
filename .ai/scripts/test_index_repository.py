#!/usr/bin/env python3
"""Tests for the dependency-free CFDX AI repository indexer."""

from __future__ import annotations

import pathlib
import sqlite3
import subprocess
import sys
import tempfile
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parent / "index_repository.py"


class IndexRepositoryTests(unittest.TestCase):
    def test_indexes_supported_sources_deterministically(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            (root / "src").mkdir()
            (root / "src" / "b.cpp").write_text("int b() { return 2; }\n", encoding="utf-8")
            (root / "src" / "a.py").write_text("x = 1\n", encoding="utf-8")
            (root / "README.md").write_text("ignored\n", encoding="utf-8")
            output = root / "index.sqlite"

            subprocess.run(
                [sys.executable, str(SCRIPT), "--root", str(root), "--output", str(output)],
                check=True,
            )
            with sqlite3.connect(output) as db:
                rows = db.execute(
                    "SELECT path, language FROM files ORDER BY path"
                ).fetchall()
                metadata = dict(db.execute("SELECT key, value FROM index_metadata"))

            self.assertEqual(rows, [("src/a.py", "python"), ("src/b.cpp", "cpp")])
            self.assertEqual(metadata["schema_version"], "1")
            self.assertEqual(metadata["repository_revision"], "unknown")

    def test_excludes_generated_directories(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            (root / "build").mkdir()
            (root / "build" / "generated.cpp").write_text("int x;\n", encoding="utf-8")
            output = root / "index.sqlite"

            subprocess.run(
                [sys.executable, str(SCRIPT), "--root", str(root), "--output", str(output)],
                check=True,
            )
            with sqlite3.connect(output) as db:
                count = db.execute("SELECT COUNT(*) FROM files").fetchone()[0]

            self.assertEqual(count, 0)


if __name__ == "__main__":
    unittest.main()
