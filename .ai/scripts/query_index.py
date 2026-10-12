#!/usr/bin/env python3
"""Query and validate a CFDX AI SQLite code index."""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import sqlite3
import subprocess
import sys


def revision(root: pathlib.Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def metadata(db: sqlite3.Connection) -> dict[str, str]:
    return dict(db.execute("SELECT key, value FROM index_metadata"))


def validate(db: sqlite3.Connection, root: pathlib.Path) -> list[str]:
    meta = metadata(db)
    errors: list[str] = []
    expected = meta.get("repository_revision")
    actual = revision(root)
    if expected and expected != "unknown" and actual != "unknown" and expected != actual:
        errors.append(f"repository revision mismatch: index={expected} current={actual}")
    for row in db.execute("SELECT path, content_hash FROM files ORDER BY path"):
        path = root / row[0]
        if not path.is_file():
            errors.append(f"indexed file missing: {row[0]}")
        elif sha256(path) != row[1]:
            errors.append(f"indexed file changed: {row[0]}")
    return errors


def symbols(db: sqlite3.Connection, pattern: str) -> list[tuple]:
    return db.execute(
        """SELECT symbols.kind, symbols.qualified_name, files.path,
                  symbols.line_start, symbols.column_start
           FROM symbols JOIN files ON files.id = symbols.file_id
           WHERE symbols.name LIKE ? OR symbols.qualified_name LIKE ?
           ORDER BY files.path, symbols.line_start, symbols.column_start""",
        (f"%{pattern}%", f"%{pattern}%"),
    ).fetchall()


def refs(db: sqlite3.Connection, pattern: str) -> list[tuple]:
    return db.execute(
        """SELECT r.kind, files.path, r.line, r.column,
                  target.name, target.qualified_name
           FROM code_references r
           JOIN symbols target ON target.id = r.target_symbol_id
           JOIN files ON files.id = r.file_id
           WHERE target.name LIKE ? OR target.qualified_name LIKE ?
           ORDER BY files.path, r.line, r.column""",
        (f"%{pattern}%", f"%{pattern}%"),
    ).fetchall()


def deps(db: sqlite3.Connection, pattern: str) -> list[tuple]:
    return db.execute(
        """SELECT files.path, dependencies.kind, dependencies.target_path
           FROM dependencies JOIN files ON files.id = dependencies.source_file_id
           WHERE files.path LIKE ? OR dependencies.target_path LIKE ?
           ORDER BY files.path, dependencies.target_path""",
        (f"%{pattern}%", f"%{pattern}%"),
    ).fetchall()


def tests(db: sqlite3.Connection, pattern: str) -> list[tuple]:
    return db.execute(
        """SELECT tests.name, tests.framework, files.path,
                  symbols.kind, symbols.qualified_name, test_links.relation
           FROM tests
           JOIN files ON files.id = tests.file_id
           LEFT JOIN test_links ON test_links.test_id = tests.id
           LEFT JOIN symbols ON symbols.id = test_links.symbol_id
           WHERE tests.name LIKE ? OR files.path LIKE ?
           ORDER BY files.path, tests.name, symbols.line_start, symbols.column_start""",
        (f"%{pattern}%", f"%{pattern}%"),
    ).fetchall()


def validations(db: sqlite3.Connection, pattern: str) -> list[tuple]:
    return db.execute(
        """SELECT validation_cases.name, validation_cases.path,
                  validation_cases.category, validation_cases.status,
                  symbols.kind, symbols.qualified_name, validation_links.relation
           FROM validation_cases
           LEFT JOIN validation_links
             ON validation_links.validation_case_id = validation_cases.id
           LEFT JOIN symbols ON symbols.id = validation_links.symbol_id
           WHERE validation_cases.name LIKE ? OR validation_cases.path LIKE ?
           ORDER BY validation_cases.path, validation_cases.name,
                    symbols.line_start, symbols.column_start""",
        (f"%{pattern}%", f"%{pattern}%"),
    ).fetchall()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=pathlib.Path, required=True)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--validate", action="store_true")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--symbol")
    group.add_argument("--references")
    group.add_argument("--dependency")
    group.add_argument("--test")
    group.add_argument("--validation")
    args = parser.parse_args()

    root = args.root.resolve()
    with sqlite3.connect(args.index) as db:
        if args.validate:
            errors = validate(db, root)
            if errors:
                for error in errors:
                    print(f"STALE: {error}", file=sys.stderr)
                return 2
            print("INDEX VALID")
            return 0

        if args.symbol:
            rows = symbols(db, args.symbol)
        elif args.references:
            rows = refs(db, args.references)
        elif args.dependency:
            rows = deps(db, args.dependency)
        elif args.test:
            rows = tests(db, args.test)
        elif args.validation:
            rows = validations(db, args.validation)
        else:
            parser.error("one query option or --validate is required")

        for row in rows:
            print("\t".join("" if value is None else str(value) for value in row))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
