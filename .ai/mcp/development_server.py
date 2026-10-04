#!/usr/bin/env python3
"""Read-only MCP adapter for the CFDX code-intelligence query contract."""
from __future__ import annotations

import argparse
import os
import pathlib
import subprocess
import sys
from typing import Annotated

from pydantic import Field
from mcp.server import MCPServer
from mcp.types import ToolAnnotations

QUERY_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "query_index.py"


def _config(root: str | None, index: str | None) -> tuple[pathlib.Path, pathlib.Path]:
    repository = pathlib.Path(root or os.environ.get("CFDX_AI_ROOT", ".")).resolve()
    database = pathlib.Path(
        index
        or os.environ.get("CFDX_AI_INDEX", ".ai/index/cfdx-ai-index.sqlite")
    ).resolve()
    return repository, database


def _run_query(root: pathlib.Path, index: pathlib.Path, *args: str) -> tuple[int, str, str]:
    if not index.is_file():
        return 2, "", f"index file does not exist: {index}"
    command = [
        sys.executable,
        str(QUERY_SCRIPT),
        "--root",
        str(root),
        "--index",
        str(index),
        *args,
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    return completed.returncode, completed.stdout, completed.stderr


def _revision(root: pathlib.Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _validate(root: pathlib.Path, index: pathlib.Path) -> dict[str, object]:
    code, stdout, stderr = _run_query(root, index, "--validate")
    if code == 0:
        return {
            "ok": True,
            "freshness": "valid",
            "repository_revision": _revision(root),
            "index": str(index),
            "errors": [],
        }
    errors = [line.removeprefix("STALE: ").strip() for line in stderr.splitlines() if line.strip()]
    if not errors and stderr.strip():
        errors = [stderr.strip()]
    return {
        "ok": False,
        "freshness": "stale" if code == 2 else "error",
        "repository_revision": _revision(root),
        "index": str(index),
        "errors": errors,
        "stdout": stdout.strip(),
    }


def _query(
    operation: str,
    pattern: str,
    root: str | None,
    index: str | None,
    option: str,
    columns: tuple[str, ...],
) -> dict[str, object]:
    repository, database = _config(root, index)
    freshness = _validate(repository, database)
    if not freshness["ok"]:
        return {
            "operation": operation,
            **freshness,
            "rows": [],
        }

    code, stdout, stderr = _run_query(repository, database, option, pattern)
    if code != 0:
        return {
            "operation": operation,
            **freshness,
            "ok": False,
            "freshness": "error",
            "errors": [stderr.strip() or "query failed"],
            "rows": [],
        }

    rows = []
    for line in stdout.splitlines():
        values = line.split("\t")
        rows.append(dict(zip(columns, values)))
    return {
        "operation": operation,
        **freshness,
        "ok": True,
        "rows": rows,
    }


def _git_ls_files(root: pathlib.Path, path: str) -> dict[str, object]:
    if path.startswith("/") or ".." in pathlib.PurePosixPath(path).parts:
        return {"ok": False, "errors": ["path must remain inside the repository"], "files": []}
    try:
        completed = subprocess.run(
            ["git", "-C", str(root), "ls-files", "--", path or "."],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        return {"ok": False, "errors": [str(exc)], "files": []}
    if completed.returncode != 0:
        return {"ok": False, "errors": [completed.stderr.strip() or "git ls-files failed"], "files": []}
    return {"ok": True, "files": completed.stdout.splitlines()}


def _git_search(root: pathlib.Path, pattern: str, path: str) -> dict[str, object]:
    if path.startswith("/") or ".." in pathlib.PurePosixPath(path).parts:
        return {"ok": False, "errors": ["path must remain inside the repository"], "matches": []}
    try:
        completed = subprocess.run(
            ["git", "-C", str(root), "grep", "-n", "-I", "--", pattern, path or "."],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        return {"ok": False, "errors": [str(exc)], "matches": []}
    if completed.returncode not in (0, 1):
        return {"ok": False, "errors": [completed.stderr.strip() or "git grep failed"], "matches": []}
    matches = []
    for line in completed.stdout.splitlines():
        parts = line.split(":", 2)
        if len(parts) == 3:
            matches.append({"path": parts[0], "line": parts[1], "text": parts[2]})
    return {"ok": True, "errors": [], "matches": matches}


def _git_read(root: pathlib.Path, path: str) -> dict[str, object]:
    if path.startswith("/") or ".." in pathlib.PurePosixPath(path).parts or not path:
        return {"ok": False, "errors": ["path must be a non-empty repository-relative tracked path"], "content": ""}
    try:
        completed = subprocess.run(
            ["git", "-C", str(root), "show", f"HEAD:{path}"],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        return {"ok": False, "errors": [str(exc)], "content": ""}
    if completed.returncode != 0:
        return {"ok": False, "errors": [completed.stderr.strip() or "tracked file could not be read"], "content": ""}
    return {"ok": True, "errors": [], "path": path, "content": completed.stdout}


def create_server(root: str | None = None, index: str | None = None) -> MCPServer:
    repository, database = _config(root, index)
    server = MCPServer(
        "CFDX Development MCP",
        version="0.1.0",
        instructions=(
            "Read-only CFDX code-intelligence server. "
            "Indexed evidence is descriptive and is never execution coverage, "
            "pass/fail, verification, or qualification evidence."
        ),
    )
    annotations = ToolAnnotations(read_only_hint=True, open_world_hint=False)

    @server.tool(
        name="index.validate",
        title="Validate code index",
        annotations=annotations,
    )
    def index_validate() -> dict[str, object]:
        """Validate the configured repository/index pair before evidence use."""
        return _validate(repository, database)

    @server.tool(
        name="code.symbol",
        title="Query symbols",
        annotations=annotations,
    )
    def code_symbol(pattern: str) -> dict[str, object]:
        """Query indexed symbols by name or qualified-name pattern."""
        return _query(
            "code.symbol",
            pattern,
            str(repository),
            str(database),
            "--symbol",
            ("kind", "qualified_name", "path", "line_start", "column_start"),
        )

    @server.tool(
        name="code.references",
        title="Query code references",
        annotations=annotations,
    )
    def code_references(pattern: str) -> dict[str, object]:
        """Query indexed references by target symbol pattern."""
        return _query(
            "code.references",
            pattern,
            str(repository),
            str(database),
            "--references",
            ("kind", "path", "line", "column", "target_name", "target_qualified_name"),
        )

    @server.tool(
        name="code.dependencies",
        title="Query dependencies",
        annotations=annotations,
    )
    def code_dependencies(pattern: str) -> dict[str, object]:
        """Query indexed dependency edges by path pattern."""
        return _query(
            "code.dependencies",
            pattern,
            str(repository),
            str(database),
            "--dependency",
            ("path", "kind", "target_path"),
        )

    @server.tool(
        name="evidence.test",
        title="Query test evidence",
        annotations=annotations,
    )
    def evidence_test(pattern: str) -> dict[str, object]:
        """Query discovered tests and explicit test-to-symbol relations."""
        return _query(
            "evidence.test",
            pattern,
            str(repository),
            str(database),
            "--test",
            ("test_name", "framework", "path", "symbol_kind", "symbol", "relation"),
        )

    @server.tool(
        name="evidence.validation",
        title="Query validation evidence",
        annotations=annotations,
    )
    def evidence_validation(pattern: str) -> dict[str, object]:
        """Query validation cases and explicit validation-to-symbol relations."""
        return _query(
            "evidence.validation",
            pattern,
            str(repository),
            str(database),
            "--validation",
            (
                "validation_name",
                "path",
                "category",
                "status",
                "symbol_kind",
                "symbol",
                "relation",
            ),
        )

    @server.tool(
        name="repository.file_structure",
        title="List tracked repository files",
        annotations=annotations,
    )
    def repository_file_structure(path: str = "") -> dict[str, object]:
        """List tracked Git files below an optional repository-relative path."""
        return _git_ls_files(repository, path)

    @server.tool(
        name="repository.status",
        title="Inspect repository status",
        annotations=annotations,
    )
    def repository_status() -> dict[str, object]:
        """Inspect the current Git branch/status without modifying repository state."""
        try:
            completed = subprocess.run(
                ["git", "-C", str(repository), "status", "--short", "--branch"],
                text=True,
                capture_output=True,
                check=False,
            )
        except OSError as exc:
            return {"ok": False, "errors": [str(exc)], "status": []}
        if completed.returncode != 0:
            return {
                "ok": False,
                "errors": [completed.stderr.strip() or "git status failed"],
                "status": [],
            }
        return {"ok": True, "errors": [], "status": completed.stdout.splitlines()}

    @server.tool(
        name="documentation.search",
        title="Search tracked documentation",
        annotations=annotations,
    )
    def documentation_search(pattern: str, path: str = "") -> dict[str, object]:
        """Search tracked repository text without modifying repository state."""
        return _git_search(repository, pattern, path)

    @server.tool(
        name="documentation.read",
        title="Read tracked documentation file",
        annotations=annotations,
    )
    def documentation_read(
        path: Annotated[str, Field(validation_alias="file_path")]
    ) -> dict[str, object]:
        """Read a tracked repository file at HEAD without modifying repository state."""
        return _git_read(repository, path)

    return server


mcp = create_server()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=None)
    parser.add_argument("--index", type=pathlib.Path, default=None)
    parser.add_argument(
        "--transport",
        choices=("stdio", "streamable-http"),
        default="stdio",
    )
    args = parser.parse_args()

    create_server(args.root, args.index).run(transport=args.transport)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
