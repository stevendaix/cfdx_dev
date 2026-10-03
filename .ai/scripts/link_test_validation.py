"""Create conservative semantic test/V&V links from explicit repository evidence.

This first backend intentionally links only evidence that is explicit in source:
- CTest/CMake test registrations to discovered test source files when the path is explicit.
- Validation manifests that explicitly name a test path.
No filename similarity, symbol proximity, or guessed coverage is used.
"""
from __future__ import annotations
import argparse
import json
import re
import sqlite3
from pathlib import Path


def files_under(root: Path):
    for p in sorted(root.rglob("*")):
        if not p.is_file() or any(part in {".git", "build", ".venv", "__pycache__"} for part in p.parts):
            continue
        yield p


def rel(p: Path, root: Path) -> str:
    return p.relative_to(root).as_posix()


def test_paths(root: Path) -> set[str]:
    out = set()
    for p in files_under(root / "tests"):
        r = rel(p, root)
        n = p.name.lower()
        if p.suffix.lower() in {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".py"} and (n.startswith("test_") or any(x in n for x in ("_test.", "_tests.", "tests."))):
            out.add(r)
    return out


def explicit_test_links(root: Path, known: set[str]):
    links = []
    # Only accept literal test source paths appearing in CMake test registration context.
    pat = re.compile(r"(?:add_test|add_executable|py_test|pytest).*", re.I)
    for p in files_under(root):
        if p.suffix.lower() not in {".cmake", ".txt"}:
            continue
        text = p.read_text(encoding="utf-8", errors="ignore")
        for line in text.splitlines():
            if not pat.search(line):
                continue
            for t in known:
                if re.search(r"(?<![A-Za-z0-9_./-])" + re.escape(t) + r"(?![A-Za-z0-9_./-])", line):
                    links.append((t, rel(p, root), "explicit-cmake-registration"))
    return sorted(set(links))


def validation_links(root: Path, known: set[str]):
    links = []
    base = root / "tests" / "validation"
    if not base.exists():
        return links
    for p in files_under(base):
        if p.suffix.lower() not in {".json", ".yaml", ".yml", ".toml", ".md"}:
            continue
        text = p.read_text(encoding="utf-8", errors="ignore")
        for t in known:
            if re.search(r"(?<![A-Za-z0-9_./-])" + re.escape(t) + r"(?![A-Za-z0-9_./-])", text):
                links.append((rel(p, root), t, "explicit-validation-reference"))
    return sorted(set(links))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--index", type=Path, required=True)
    args = ap.parse_args()
    root = args.root.resolve(); index = args.index.resolve()
    known = test_paths(root)
    conn = sqlite3.connect(index)
    # Keep this backend deliberately non-authoritative: test_links only stores explicit evidence.
    test_id = {}
    for path in known:
        row = conn.execute("SELECT id FROM tests WHERE file_id IN (SELECT id FROM files WHERE path=?)", (path,)).fetchone()
        if row: test_id[path] = row[0]
    validation_id = {r["path"]: r["id"] for r in conn.execute("SELECT id,path FROM validation_cases").fetchall()}

    # A test source explicitly defines the symbols in that same source file.
    # This is traceability, not execution/coverage proof.
    for path, tid in test_id.items():
        rows = conn.execute("SELECT s.id FROM symbols s JOIN files f ON f.id=s.file_id WHERE f.path=?", (path,)).fetchall()
        for (sid,) in rows:
            conn.execute("INSERT OR IGNORE INTO test_links(test_id, symbol_id, relation) VALUES (?, ?, ?)", (tid, sid, "defines-symbol"))

    for vpath, tpath, evidence in validation_links(root, known):
        vid = validation_id.get(vpath)
        if vid and tpath in test_id:
            rows = conn.execute("SELECT s.id FROM symbols s JOIN files f ON f.id=s.file_id WHERE f.path=?", (tpath,)).fetchall()
            for (sid,) in rows:
                conn.execute("INSERT OR IGNORE INTO validation_links(validation_case_id, symbol_id, relation) VALUES (?, ?, ?)", (vid, sid, "explicit-test-reference"))
    conn.execute("INSERT OR REPLACE INTO index_metadata(key,value) VALUES(?,?)", ("semantic_linkage", "explicit-evidence-only-v1"))
    conn.commit(); conn.close()

if __name__ == "__main__":
    main()
