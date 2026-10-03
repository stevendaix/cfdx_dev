from __future__ import annotations
import sqlite3, tempfile
from pathlib import Path
import subprocess, sys

def test_explicit_only():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root/"tests/validation").mkdir(parents=True)
        (root/"tests/unit").mkdir(parents=True)
        (root/"tests/unit/test_alpha.cpp").write_text("int alpha_test() { return 0; }\n")
        (root/"tests/validation/case.yml").write_text("test: tests/unit/test_alpha.cpp\n")
        (root/"CMakeLists.txt").write_text("add_test(NAME alpha COMMAND tests/unit/test_alpha.cpp)\n")
        db = root/"index.sqlite"
        c=sqlite3.connect(db)
        c.executescript("""CREATE TABLE files(id INTEGER PRIMARY KEY,path TEXT UNIQUE); CREATE TABLE tests(id INTEGER PRIMARY KEY,file_id INTEGER); CREATE TABLE symbols(id INTEGER PRIMARY KEY,file_id INTEGER); CREATE TABLE validation_cases(id INTEGER PRIMARY KEY,path TEXT UNIQUE); CREATE TABLE test_links(test_id INTEGER,symbol_id INTEGER,relation TEXT, UNIQUE(test_id,symbol_id,relation)); CREATE TABLE validation_links(validation_case_id INTEGER,symbol_id INTEGER,relation TEXT, UNIQUE(validation_case_id,symbol_id,relation)); CREATE TABLE index_metadata(key TEXT PRIMARY KEY,value TEXT); INSERT INTO files(id,path) VALUES(1,'tests/unit/test_alpha.cpp'); INSERT INTO tests(id,file_id) VALUES(1,1); INSERT INTO symbols(id,file_id) VALUES(1,1); INSERT INTO validation_cases(id,path) VALUES(1,'tests/validation/case.yml');""")
        c.commit(); c.close()
        subprocess.run([sys.executable, ".ai/scripts/link_test_validation.py", "--root", str(root), "--index", str(db)], check=True)
        c=sqlite3.connect(db)
        assert c.execute("SELECT relation FROM test_links").fetchall() == [("defines-symbol",)]
        assert c.execute("SELECT relation FROM validation_links").fetchall() == [("explicit-test-reference",)]
        c.close()

if __name__ == "__main__":
    test_explicit_only()
