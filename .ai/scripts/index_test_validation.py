#!/usr/bin/env python3
"""Discover repository test sources and validation cases in the AI index.

This inventory is deliberately conservative: it records discoverable repository
artifacts and provenance, but does not infer semantic test-to-symbol or
validation-to-symbol coverage.
"""
from __future__ import annotations

import argparse
import pathlib
import sqlite3

SOURCE_SUFFIXES = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".py"}
TEST_ROOTS = ("tests",)
TEST_NAME_MARKERS = ("test_", "_test.", "_tests.", "tests.")
VALIDATION_ROOT = pathlib.Path("tests/validation")


def is_test_source(path: pathlib.Path) -> bool:
    if path.suffix.lower() not in SOURCE_SUFFIXES:
        return False
    name = path.name.lower()
    return name.startswith("test_") or any(marker in name for marker in TEST_NAME_MARKERS)


def discover_tests(root: pathlib.Path) -> list[tuple[str, str]]:
    found: list[tuple[str, str]] = []
    for root_name in TEST_ROOTS:
        base = root / root_name
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and is_test_source(path):
                found.append((path.relative_to(root).as_posix(), "source-discovery"))
    return found


def discover_validation_cases(root: pathlib.Path) -> list[tuple[str, str, str]]:
    base = root / VALIDATION_ROOT
    if not base.is_dir():
        return []
    cases: list[tuple[str, str, str]] = []
    for path in sorted(base.iterdir()):
        if path.name.startswith("."):
            continue
        if path.is_dir():
            cases.append((path.relative_to(root).as_posix(), path.name, "directory-discovery"))
        elif path.is_file() and path.suffix.lower() in {".json", ".yaml", ".yml", ".toml", ".md"}:
            cases.append((path.relative_to(root).as_posix(), path.stem, "manifest-discovery"))
    return cases


def index_test_validation(root: pathlib.Path, db_path: pathlib.Path) -> tuple[int, int]:
    root = root.resolve()
    with sqlite3.connect(db_path) as db:
        db.execute("PRAGMA foreign_keys = ON")
        db.execute(
            "INSERT INTO index_metadata(key, value) VALUES (?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("test_validation_discovery", "conservative-source-and-validation-inventory-v1"),
        )

        test_count = 0
        for relpath, provenance in discover_tests(root):
            row = db.execute("SELECT id FROM files WHERE path = ?", (relpath,)).fetchone()
            if row is None:
                continue
            db.execute(
                "INSERT INTO tests(file_id, name, framework, line_start, line_end) VALUES (?, ?, ?, ?, ?)",
                (row[0], pathlib.Path(relpath).stem, provenance, None, None),
            )
            test_count += 1

        validation_count = 0
        for relpath, name, provenance in discover_validation_cases(root):
            db.execute(
                "INSERT OR IGNORE INTO validation_cases(path, name, category, status, source_revision) "
                "VALUES (?, ?, ?, ?, ?)",
                (relpath, name, provenance, "discovered", None),
            )
            validation_count += 1

    return test_count, validation_count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--index", type=pathlib.Path, required=True)
    args = parser.parse_args()
    tests, validations = index_test_validation(args.root, args.index)
    print(f"DISCOVERED tests={tests} validation_cases={validations}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
