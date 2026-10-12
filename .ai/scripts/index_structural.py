#!/usr/bin/env python3
"""Populate structural symbols/dependencies in a CFDX AI SQLite index.

Python uses the standard-library AST. C/C++ currently uses a conservative
declaration scanner and is explicitly marked as heuristic; it never claims
complete semantic C++ coverage.
"""
from __future__ import annotations

import ast
import pathlib
import re
import sqlite3
from dataclasses import dataclass

CPP_DECL = re.compile(
    r"(?m)^\s*(?:(?:template\s*<[^;{}]+>\s*)?"
    r"(?:inline\s+|static\s+|virtual\s+|constexpr\s+|const\s+|extern\s+)*"
    r"(?:[A-Za-z_][\w:<>]*\s+)+)"
    r"([A-Za-z_]\w*)\s*\([^;{}]*\)\s*(?:const\s*)?(?:\{|;)"
)
CPP_INCLUDE = re.compile(r'(?m)^\s*#\s*include\s*[<"]([^>"]+)[>"]')


@dataclass(frozen=True)
class Symbol:
    kind: str
    qualified_name: str
    name: str
    line_start: int
    column_start: int | None
    line_end: int | None
    column_end: int | None


def _python_symbols(source: str) -> tuple[list[Symbol], list[tuple[str, str]]]:
    tree = ast.parse(source)
    symbols: list[Symbol] = []
    dependencies: list[tuple[str, str]] = []

    def visit(node: ast.AST, prefix: str = "") -> None:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            q = f"{prefix}.{node.name}" if prefix else node.name
            symbols.append(Symbol(
                "function", q, node.name, node.lineno,
                node.col_offset, getattr(node, "end_lineno", None),
                getattr(node, "end_col_offset", None)))
            for child in node.body:
                visit(child, q)
            return
        if isinstance(node, ast.ClassDef):
            q = f"{prefix}.{node.name}" if prefix else node.name
            symbols.append(Symbol(
                "class", q, node.name, node.lineno,
                node.col_offset, getattr(node, "end_lineno", None),
                getattr(node, "end_col_offset", None)))
            for child in node.body:
                visit(child, q)
            return
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    dependencies.append(("python-import", alias.name))
            else:
                module = node.module or ""
                dependencies.append(("python-import", module))
        for child in ast.iter_child_nodes(node):
            visit(child, prefix)

    visit(tree)
    return symbols, dependencies


def _cpp_symbols(source: str) -> tuple[list[Symbol], list[tuple[str, str]]]:
    symbols = [
        Symbol("function", name, name, match.start(1) and source.count("\n", 0, match.start(1)) + 1,
               match.start(1), None, None)
        for match in CPP_DECL.finditer(source)
        for name in [match.group(1)]
    ]
    dependencies = [("cpp-include", match.group(1)) for match in CPP_INCLUDE.finditer(source)]
    return symbols, dependencies


def index_structural(root: pathlib.Path, db_path: pathlib.Path) -> None:
    root = root.resolve()
    with sqlite3.connect(db_path) as db:
        db.execute("PRAGMA foreign_keys = ON")
        file_rows = db.execute(
            "SELECT id, path, language FROM files ORDER BY path"
        ).fetchall()
        for file_id, relpath, language in file_rows:
            source = (root / relpath).read_text(encoding="utf-8", errors="replace")
            try:
                if language == "python":
                    symbols, dependencies = _python_symbols(source)
                    parser = "python-ast"
                elif language in {"c", "cpp", "cpp-header"}:
                    symbols, dependencies = _cpp_symbols(source)
                    parser = "cpp-declaration-scanner"
                else:
                    continue
            except (SyntaxError, UnicodeError) as exc:
                db.execute(
                    "INSERT INTO index_metadata(key, value) VALUES (?, ?) "
                    "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                    (f"parse_error:{relpath}", f"{type(exc).__name__}: {exc}"),
                )
                continue

            db.execute(
                "INSERT INTO index_metadata(key, value) VALUES (?, ?) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (f"parser:{relpath}", parser),
            )
            for symbol in symbols:
                db.execute(
                    """INSERT INTO symbols(
                        file_id, kind, qualified_name, name, line_start,
                        column_start, line_end, column_end)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                    (file_id, symbol.kind, symbol.qualified_name, symbol.name,
                     symbol.line_start, symbol.column_start, symbol.line_end,
                     symbol.column_end),
                )
            for kind, target in dependencies:
                db.execute(
                    "INSERT INTO dependencies(source_file_id, target_path, kind) "
                    "VALUES (?, ?, ?)",
                    (file_id, target, kind),
                )


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--index", type=pathlib.Path, required=True)
    args = parser.parse_args()
    index_structural(args.root, args.index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
