#!/usr/bin/env python3
"""Index C/C++ semantics from a CMake compilation database using libclang.

The backend is optional by design. It never pretends that the heuristic
scanner is semantic C++ and records parser provenance/diagnostics explicitly.
"""
from __future__ import annotations

import json
import pathlib
import shlex
import sqlite3
from dataclasses import dataclass


@dataclass(frozen=True)
class CompileCommand:
    directory: pathlib.Path
    file: pathlib.Path
    arguments: tuple[str, ...]


def load_compile_commands(path: pathlib.Path) -> list[CompileCommand]:
    entries = json.loads(path.read_text(encoding="utf-8"))
    commands: list[CompileCommand] = []
    for entry in entries:
        directory = pathlib.Path(entry["directory"]).resolve()
        source = pathlib.Path(entry["file"])
        if not source.is_absolute():
            source = directory / source
        raw = entry.get("arguments")
        arguments = tuple(raw) if raw else tuple(shlex.split(entry["command"]))
        commands.append(CompileCommand(directory, source.resolve(), arguments))
    return commands


def _load_clang():
    try:
        from clang import cindex
    except ImportError as exc:
        raise RuntimeError("libclang Python bindings are unavailable") from exc
    return cindex


def _clang_arguments(command: CompileCommand) -> list[str]:
    args = list(command.arguments)
    if args and pathlib.Path(args[0]).name in {"cc", "gcc", "clang", "c++", "g++", "clang++"}:
        args = args[1:]
    filtered: list[str] = []
    skip = False
    for arg in args:
        if skip:
            skip = False
            continue
        if arg in {"-c", "-o"}:
            skip = arg == "-o"
            continue
        if arg == str(command.file) or pathlib.Path(arg).resolve() == command.file:
            continue
        filtered.append(arg)
    return filtered


def index_clang(root: pathlib.Path, db_path: pathlib.Path, compile_commands: pathlib.Path) -> int:
    cindex = _load_clang()
    commands = load_compile_commands(compile_commands)
    root = root.resolve()
    indexed = 0

    with sqlite3.connect(db_path) as db:
        db.execute("PRAGMA foreign_keys = ON")
        for command in commands:
            try:
                relpath = command.file.relative_to(root).as_posix()
            except ValueError:
                continue
            row = db.execute("SELECT id FROM files WHERE path = ?", (relpath,)).fetchone()
            if row is None:
                continue

            try:
                index = cindex.Index.create()
                tu = index.parse(
                    str(command.file),
                    args=_clang_arguments(command),
                    options=0,
                )
            except Exception as exc:
                db.execute(
                    "INSERT INTO index_metadata(key, value) VALUES (?, ?) "
                    "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                    (f"clang_error:{relpath}", f"{type(exc).__name__}: {exc}"),
                )
                continue

            db.execute(
                "INSERT INTO index_metadata(key, value) VALUES (?, ?) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (f"parser:{relpath}", "clang-libclang"),
            )
            db.execute(
                "INSERT INTO index_metadata(key, value) VALUES (?, ?) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (f"coverage:{relpath}", "semantic-cpp"),
            )

            file_id = row[0]
            symbol_ids: dict[tuple[str, int, int], int] = {}

            def visit(cursor, owner: int | None = None) -> None:
                loc = cursor.location
                if not loc.file or pathlib.Path(loc.file.name).resolve() != command.file:
                    for child in cursor.get_children():
                        visit(child, owner)
                    return

                spelling = cursor.spelling or cursor.displayname
                kind_name = str(cursor.kind).split(".")[-1].lower()
                symbol_kind = None
                if cursor.kind in {
                    cindex.CursorKind.NAMESPACE,
                    cindex.CursorKind.CLASS_DECL,
                    cindex.CursorKind.STRUCT_DECL,
                    cindex.CursorKind.ENUM_DECL,
                    cindex.CursorKind.FUNCTION_DECL,
                    cindex.CursorKind.FUNCTION_TEMPLATE,
                    cindex.CursorKind.CXX_METHOD,
                    cindex.CursorKind.CONSTRUCTOR,
                    cindex.CursorKind.DESTRUCTOR,
                    cindex.CursorKind.VAR_DECL,
                } and spelling:
                    symbol_kind = kind_name
                    qname = cursor.displayname or spelling
                    result = db.execute(
                        """INSERT INTO symbols(
                           file_id, kind, qualified_name, name, line_start,
                           column_start, line_end, column_end)
                           VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                        (file_id, symbol_kind, qname, spelling, loc.line,
                         loc.column, None, None),
                    )
                    symbol_ids[(spelling, loc.line, loc.column)] = result.lastrowid
                    owner = result.lastrowid

                if cursor.kind == cindex.CursorKind.INCLUSION_DIRECTIVE:
                    target = cursor.displayname or spelling
                    if target:
                        db.execute(
                            "INSERT INTO dependencies(source_file_id, target_path, kind) VALUES (?, ?, ?)",
                            (file_id, target, "clang-include"),
                        )

                if cursor.kind in {
                    cindex.CursorKind.CALL_EXPR,
                    cindex.CursorKind.DECL_REF_EXPR,
                    cindex.CursorKind.MEMBER_REF_EXPR,
                }:
                    referenced = cursor.referenced
                    if referenced is not None and referenced.location.file:
                        ref_file = pathlib.Path(referenced.location.file.name).resolve()
                        if ref_file == command.file:
                            target = symbol_ids.get(
                                (referenced.spelling, referenced.location.line, referenced.location.column)
                            )
                            if target is not None:
                                db.execute(
                                    """INSERT INTO references(
                                       source_symbol_id, target_symbol_id, file_id,
                                       kind, line, column)
                                       VALUES (?, ?, ?, ?, ?, ?)""",
                                    (owner, target, kind_name, file_id, loc.line, loc.column),
                                )

                for child in cursor.get_children():
                    visit(child, owner)

            visit(tu.cursor)
            diagnostics = [str(d) for d in tu.diagnostics if d.severity >= 3]
            if diagnostics:
                db.execute(
                    "INSERT INTO index_metadata(key, value) VALUES (?, ?) "
                    "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                    (f"clang_diagnostics:{relpath}", "\n".join(diagnostics)),
                )
            indexed += 1

    return indexed


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--index", type=pathlib.Path, required=True)
    parser.add_argument("--compile-commands", type=pathlib.Path, required=True)
    args = parser.parse_args()
    try:
        index_clang(args.root, args.index, args.compile_commands)
    except RuntimeError as exc:
        print(f"CFDX clang indexer unavailable: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
