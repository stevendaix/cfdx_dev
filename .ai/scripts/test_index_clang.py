#!/usr/bin/env python3
"""Tests for the optional compile_commands/libclang backend."""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sqlite3
import tempfile
import unittest


SCRIPT = pathlib.Path(__file__).with_name("index_clang.py")
SPEC = importlib.util.spec_from_file_location("index_clang", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ClangIndexerTests(unittest.TestCase):
    def test_compile_commands_arguments_and_relative_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            source = root / "src" / "example.cpp"
            source.parent.mkdir()
            source.write_text("int f() { return 1; }\n", encoding="utf-8")
            database = root / "compile_commands.json"
            database.write_text(json.dumps([{
                "directory": str(root),
                "file": "src/example.cpp",
                "arguments": ["clang++", "-std=c++17", "-c", "src/example.cpp", "-o", "example.o"],
            }]), encoding="utf-8")
            commands = MODULE.load_compile_commands(database)
            self.assertEqual(commands[0].file, source.resolve())
            self.assertEqual(MODULE._clang_arguments(commands[0]), ["-std=c++17"])

    def test_command_form_is_supported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            source = root / "x.cpp"
            database = root / "compile_commands.json"
            database.write_text(json.dumps([{
                "directory": str(root),
                "file": "x.cpp",
                "command": "clang++ -std=c++17 -c x.cpp -o x.o",
            }]), encoding="utf-8")
            commands = MODULE.load_compile_commands(database)
            self.assertEqual(commands[0].arguments[0], "clang++")

    def test_missing_libclang_is_explicit(self):
        try:
            MODULE._load_clang()
        except RuntimeError as exc:
            self.assertIn("libclang Python bindings are unavailable", str(exc))


if __name__ == "__main__":
    unittest.main()
