#!/usr/bin/env python3
"""End-to-end tests for the dependency-free CFDX AI index pipeline."""

from __future__ import annotations

import pathlib
import sqlite3
import subprocess
import sys
import tempfile
import unittest


SCRIPTS = pathlib.Path(__file__).resolve().parent
INDEXER = SCRIPTS / "index_repository.py"
STRUCTURAL = SCRIPTS / "index_structural.py"
QUERY = SCRIPTS / "query_index.py"


class IndexPipelineTests(unittest.TestCase):
    def test_repository_structural_query_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            source = root / "solver.py"
            source.write_text(
                "def solve(value):\n    return value + 1\n",
                encoding="utf-8",
            )
            db_path = root / "index.sqlite"

            subprocess.run(
                [
                    sys.executable,
                    str(INDEXER),
                    "--root",
                    str(root),
                    "--output",
                    str(db_path),
                ],
                check=True,
            )
            subprocess.run(
                [
                    sys.executable,
                    str(STRUCTURAL),
                    "--root",
                    str(root),
                    "--index",
                    str(db_path),
                ],
                check=True,
            )
            validation = subprocess.run(
                [
                    sys.executable,
                    str(QUERY),
                    "--root",
                    str(root),
                    "--index",
                    str(db_path),
                    "--validate",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(validation.stdout.strip(), "INDEX VALID")

            query = subprocess.run(
                [
                    sys.executable,
                    str(QUERY),
                    "--root",
                    str(root),
                    "--index",
                    str(db_path),
                    "--symbol",
                    "solve",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("solver.py", query.stdout)
            self.assertIn("solve", query.stdout)

            with sqlite3.connect(db_path) as db:
                parser = db.execute(
                    "SELECT value FROM index_metadata WHERE key='parser:solver.py'"
                ).fetchone()
            self.assertEqual(parser, ("python-ast",))


if __name__ == "__main__":
    unittest.main()
