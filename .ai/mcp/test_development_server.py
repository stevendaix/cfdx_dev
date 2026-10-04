from __future__ import annotations

import asyncio
import pathlib
import sqlite3
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from mcp import Client

from development_server import create_server


def make_index(root: pathlib.Path, index: pathlib.Path) -> None:
    source = root / "src" / "demo.cpp"
    source.parent.mkdir(parents=True)
    source.write_text("int demo() { return 1; }\n", encoding="utf-8")
    subprocess.run(["git", "init", "--quiet", str(root)], check=True)
    subprocess.run(["git", "-C", str(root), "config", "user.email", "cfdx-ai@example.invalid"], check=True)
    subprocess.run(["git", "-C", str(root), "config", "user.name", "CFDX AI Test"], check=True)
    subprocess.run(["git", "-C", str(root), "add", "src/demo.cpp"], check=True)
    subprocess.run(["git", "-C", str(root), "commit", "--quiet", "-m", "fixture"], check=True)
    with sqlite3.connect(index) as db:
        db.executescript(
            """
            CREATE TABLE index_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE files (id INTEGER PRIMARY KEY, path TEXT NOT NULL, content_hash TEXT NOT NULL);
            CREATE TABLE symbols (id INTEGER PRIMARY KEY, file_id INTEGER NOT NULL, kind TEXT NOT NULL, name TEXT NOT NULL, qualified_name TEXT NOT NULL, line_start INTEGER NOT NULL, column_start INTEGER NOT NULL);
            CREATE TABLE code_references (id INTEGER PRIMARY KEY, file_id INTEGER NOT NULL, target_symbol_id INTEGER, kind TEXT NOT NULL, line INTEGER NOT NULL, column INTEGER NOT NULL);
            CREATE TABLE dependencies (id INTEGER PRIMARY KEY, source_file_id INTEGER NOT NULL, kind TEXT NOT NULL, target_path TEXT NOT NULL);
            CREATE TABLE tests (id INTEGER PRIMARY KEY, file_id INTEGER NOT NULL, name TEXT NOT NULL, framework TEXT NOT NULL);
            CREATE TABLE validation_cases (id INTEGER PRIMARY KEY, name TEXT NOT NULL, path TEXT NOT NULL, category TEXT NOT NULL, status TEXT NOT NULL);
            CREATE TABLE test_links (test_id INTEGER NOT NULL, symbol_id INTEGER NOT NULL, relation TEXT NOT NULL, PRIMARY KEY (test_id, symbol_id));
            CREATE TABLE validation_links (validation_case_id INTEGER NOT NULL, symbol_id INTEGER NOT NULL, relation TEXT NOT NULL, PRIMARY KEY (validation_case_id, symbol_id));
            """
        )
        db.execute("INSERT INTO index_metadata VALUES (?, ?)", ("repository_revision", "unknown"))
        db.execute("INSERT INTO files VALUES (?, ?, ?)", (1, "src/demo.cpp", "invalid-hash"))
        db.execute("INSERT INTO symbols VALUES (?, ?, ?, ?, ?, ?, ?)", (1, 1, "function", "demo", "demo", 1, 5))
        db.execute("INSERT INTO tests VALUES (?, ?, ?, ?)", (1, 1, "demo_test", "native"))
        db.execute("INSERT INTO validation_cases VALUES (?, ?, ?, ?, ?)", (1, "demo_validation", "tests/validation/demo", "unit", "planned"))


async def exercise() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        root = pathlib.Path(directory)
        index = root / "index.sqlite"
        make_index(root, index)
        client = Client(create_server(str(root), str(index)))
        async with client:
            listed = await client.list_tools()
            names = {tool.name for tool in listed.tools}
            expected_names = {
                "index.validate", "code.symbol", "code.references", "code.dependencies",
                "evidence.test", "evidence.validation", "repository.file_structure",
                "repository.status", "documentation.search", "documentation.read",
            }
            pattern_names = expected_names - {"index.validate", "repository.file_structure", "repository.status"}
            assert names == expected_names
            for tool in listed.tools:
                annotations = tool.model_dump(by_alias=True).get("annotations", {})
                assert annotations["readOnlyHint"] is True
                assert annotations["openWorldHint"] is False

            structure = await client.call_tool("repository.file_structure", {"path": "src"})
            assert structure.is_error is False
            assert structure.structured_content["ok"] is True
            assert structure.structured_content["files"] == ["src/demo.cpp"]

            status = await client.call_tool("repository.status", {})
            assert status.is_error is False
            assert status.structured_content["ok"] is True
            assert any(line.startswith("## ") for line in status.structured_content["status"])

            doc_search = await client.call_tool("documentation.search", {"pattern": "demo", "path": "src"})
            assert doc_search.is_error is False
            assert doc_search.structured_content["ok"] is True
            assert doc_search.structured_content["matches"][0]["path"] == "src/demo.cpp"

            doc_read = await client.call_tool("documentation.read", {"file_path": "src/demo.cpp"})
            assert doc_read.is_error is False
            assert doc_read.structured_content["ok"] is True
            assert "int demo()" in doc_read.structured_content["content"]

            doc_escape = await client.call_tool("documentation.read", {"file_path": "../outside"})
            assert doc_escape.is_error is False
            assert doc_escape.structured_content["ok"] is False

            doc_search_escape = await client.call_tool("documentation.search", {"pattern": "demo", "path": "../outside"})
            assert doc_search_escape.is_error is False
            assert doc_search_escape.structured_content["ok"] is False

            doc_absolute = await client.call_tool("documentation.read", {"file_path": "/etc/passwd"})
            assert doc_absolute.is_error is False
            assert doc_absolute.structured_content["ok"] is False

            escaped = await client.call_tool("repository.file_structure", {"path": "../outside"})
            assert escaped.is_error is False
            assert escaped.structured_content["ok"] is False
            assert escaped.structured_content["files"] == []

            stale = await client.call_tool("index.validate", {})
            assert stale.is_error is False
            assert stale.structured_content["freshness"] == "stale"

            for name in sorted(pattern_names):
                result = await client.call_tool(name, {"pattern": "demo"})
                assert result.is_error is False
                assert result.structured_content["freshness"] == "stale"
                assert result.structured_content["rows"] == []


if __name__ == "__main__":
    asyncio.run(exercise())
