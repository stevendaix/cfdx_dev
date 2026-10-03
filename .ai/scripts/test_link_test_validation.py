from __future__ import annotations
import sqlite3, tempfile
from pathlib import Path
import subprocess, sys


def test_explicit_only():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root/"tests/validation").mkdir(parents=True)
        (root/"tests/unit").mkdir(parents=True)
        (root/"tests/unit/test_alpha.cpp").write_text("int main(){}\n")
        (root/"tests/validation/case.yml").write_text("test: tests/unit/test_alpha.cpp\n")
        (root/"CMakeLists.txt").write_text("add_test(NAME alpha COMMAND tests/unit/test_alpha.cpp)\n")
        db = root/"index.sqlite"
        c=sqlite3.connect(db); c.executescript("""CREATE TABLE files(id INTEGER PRIMARY KEY,path TEXT UNIQUE); CREATE TABLE tests(id INTEGER PRIMARY KEY,file_id INTEGER); CREATE TABLE validation_cases(id INTEGER PRIMARY KEY,path TEXT UNIQUE); CREATE TABLE test_links(test_id INTEGER,symbol_id INTEGER,link_type TEXT, UNIQUE(test_id,symbol_id,link_type)); CREATE TABLE index_metadata(key TEXT PRIMARY KEY,value TEXT); INSERT INTO files(path) VALUES('tests/unit/test_alpha.cpp'); INSERT INTO tests(file_id) VALUES(1); INSERT INTO validation_cases(path) VALUES('tests/validation/case.yml');"""); c.commit(); c.close()
        subprocess.run([sys.executable, ".ai/scripts/link_test_validation.py", "--root", str(root), "--index", str(db)], check=True)
        c=sqlite3.connect(db); rows=c.execute("SELECT link_type FROM test_links").fetchall(); assert rows and all("explicit" in r[0] for r in rows); c.close()

if __name__ == "__main__": test_explicit_only()
