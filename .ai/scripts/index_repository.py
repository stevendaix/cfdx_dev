#!/usr/bin/env python3
"""Build a deterministic, dependency-free CFDX AI repository index."""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import sqlite3
import subprocess
from datetime import datetime, timezone

SCHEMA = pathlib.Path(__file__).resolve().parents[1] / "index" / "schema.sql"

LANGUAGES = {
    ".c": "c",
    ".cc": "cpp",
    ".cpp": "cpp",
    ".cxx": "cpp",
    ".h": "cpp-header",
    ".hh": "cpp-header",
    ".hpp": "cpp-header",
    ".hxx": "cpp-header",
    ".py": "python",
    ".cmake": "cmake",
}

EXCLUDED = {".git", "build", ".venv", "__pycache__"}


def repository_revision(root: pathlib.Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def iter_source_files(root: pathlib.Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative_parts = path.relative_to(root).parts
        if any(part in EXCLUDED for part in relative_parts):
            continue
        language = LANGUAGES.get(path.suffix.lower())
        if language is None:
            continue
        yield path, language


def content_hash(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_index(root: pathlib.Path, output: pathlib.Path) -> None:
    root = root.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()

    revision = repository_revision(root)
    now = datetime.now(timezone.utc).isoformat()
    with sqlite3.connect(output) as db:
        db.executescript(SCHEMA.read_text(encoding="utf-8"))
        db.executemany(
            "INSERT INTO index_metadata(key, value) VALUES (?, ?)",
            [
                ("schema_version", "1"),
                ("repository_revision", revision),
                ("repository_root", str(root)),
            ],
        )
        rows = []
        for path, language in iter_source_files(root):
            relative = path.relative_to(root).as_posix()
            rows.append(
                (
                    relative,
                    language,
                    content_hash(path),
                    path.stat().st_size,
                    revision,
                    now,
                )
            )
        db.executemany(
            """
            INSERT INTO files(path, language, content_hash, size_bytes, revision, indexed_at)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            rows,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument(
        "--output", type=pathlib.Path, required=True, help="Output SQLite index"
    )
    args = parser.parse_args()
    build_index(args.root, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
