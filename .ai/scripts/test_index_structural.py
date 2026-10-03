from __future__ import annotations

import pathlib
import sqlite3
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from index_repository import build_index
from index_structural import index_structural


def test_python_ast_symbols_and_imports() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        (root / "pkg").mkdir()
        (root / "pkg" / "mod.py").write_text(
            "import json\n\nclass Solver:\n    def run(self):\n        return 1\n",
            encoding="utf-8",
        )
        db_path = root / "index.sqlite"
        build_index(root, db_path)
        index_structural(root, db_path)
        with sqlite3.connect(db_path) as db:
            symbols = db.execute(
                "SELECT kind, qualified_name FROM symbols ORDER BY qualified_name"
            ).fetchall()
            deps = db.execute(
                "SELECT kind, target_path FROM dependencies"
            ).fetchall()
        assert symbols == [("class", "Solver"), ("function", "Solver.run")]
        assert deps == [("python-import", "json")]


def test_cpp_scanner_is_explicitly_heuristic() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        (root / "main.cpp").write_text(
            '#include "solver.hpp"\nint solve(int x) { return x; }\n',
            encoding="utf-8",
        )
        db_path = root / "index.sqlite"
        build_index(root, db_path)
        index_structural(root, db_path)
        with sqlite3.connect(db_path) as db:
            assert db.execute(
                "SELECT qualified_name FROM symbols"
            ).fetchall() == [("solve",)]
            assert db.execute(
                "SELECT target_path FROM dependencies"
            ).fetchall() == [("solver.hpp",)]
            assert db.execute(
                "SELECT value FROM index_metadata WHERE key='parser:main.cpp'"
            ).fetchone() == ("cpp-declaration-scanner",)
