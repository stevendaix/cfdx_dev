"""Reproducibility manifests must never invent revision/compiler provenance."""
from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "scripts" / "reproducibility_manifest.py"
spec = importlib.util.spec_from_file_location("reproducibility_manifest", SRC)
assert spec is not None and spec.loader is not None
manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(manifest)


def test_unknown_git_archive_is_not_falsely_dirty(tmp_path: Path) -> None:
    assert manifest.git_metadata(tmp_path) == {"commit": None, "dirty": None}


def test_compiler_is_from_actual_cmake_cache(tmp_path: Path) -> None:
    cxx = shutil.which("c++")
    assert cxx
    (tmp_path / "CMakeCache.txt").write_text(
        "//C++ compiler\nCMAKE_CXX_COMPILER:FILEPATH=" + cxx + "\n", encoding="utf-8"
    )
    compiler = manifest.configured_compiler(tmp_path)
    assert compiler["source"] == "CMakeCache.txt"
    assert compiler["cxx"] == cxx
    assert "version" in compiler


def test_missing_build_cache_is_rejected(tmp_path: Path) -> None:
    try:
        manifest.configured_compiler(tmp_path)
    except ValueError as exc:
        assert "cache not found" in str(exc)
    else:
        raise AssertionError("missing cache accepted")


def test_manifest_reports_hashed_case(tmp_path: Path) -> None:
    case = tmp_path / "small.case"
    case.write_text("laminar test case\n", encoding="utf-8")
    output = tmp_path / "manifest.json"
    p = subprocess.run([sys.executable, str(SRC), "--case", str(case),
                        "--output", str(output)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["case"]["sha256"] == manifest.file_sha256(case)
    assert result["compiler"]["source"] in ("PATH", "CXX")
