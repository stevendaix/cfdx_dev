"""Emit a machine-readable CFDX reproducibility manifest (metadata, not a PASS gate)."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path


def git_value(root: Path, *args: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def git_metadata(root: Path) -> dict[str, str | bool | None]:
    # An exported source ZIP is not necessarily a Git checkout. Unknown must
    # never be mislabeled as dirty=True or as a verifiable revision.
    commit = git_value(root, "rev-parse", "HEAD")
    if not isinstance(commit, str) or not re.fullmatch(r"[0-9a-f]{40}", commit):
        return {"commit": None, "dirty": None}
    status = git_value(root, "status", "--porcelain", "--untracked-files=normal")
    return {"commit": commit, "dirty": None if status is None else bool(status)}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tool_version(command: str, *args: str) -> str:
    executable = shutil.which(command)
    if not executable:
        return "unavailable"
    try:
        return (
            subprocess.check_output(
                [executable, *args], text=True, stderr=subprocess.STDOUT
            )
            .strip()
            .splitlines()[0]
        )
    except (OSError, subprocess.CalledProcessError, IndexError):
        return "unknown"


def configured_compiler(build_dir: Path | None = None) -> dict[str, str | None]:
    compiler = None
    source = "PATH"
    if build_dir is not None:
        cache = build_dir / "CMakeCache.txt"
        if not cache.is_file():
            raise ValueError(f"CMake cache not found: {cache}")
        for line in cache.read_text(encoding="utf-8").splitlines():
            if line.startswith(
                ("CMAKE_CXX_COMPILER:FILEPATH=", "CMAKE_CXX_COMPILER:STRING=")
            ):
                compiler = line.split("=", 1)[1].strip()
                source = "CMakeCache.txt"
                break
        if compiler is None:
            raise ValueError(f"CMAKE_CXX_COMPILER absent from {cache}")
    if compiler is None:
        compiler = os.environ.get("CXX") or shutil.which("c++")
        source = "CXX" if os.environ.get("CXX") else "PATH"
    return {
        "cxx": compiler,
        "version": tool_version(compiler, "--version") if compiler else "unavailable",
        "source": source,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path)
    parser.add_argument(
        "--build-dir", type=Path, help="read the real CMake C++ compiler"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    manifest = {
        "generated_at": (
            datetime.datetime.now(datetime.timezone.utc)
            .isoformat()
            .replace("+00:00", "Z")
        ),
        "cfdx_version": "0.7",
        "schema_version": 1,
        "git": git_metadata(root),
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "compiler": configured_compiler(args.build_dir),
        "toolchain": {
            "cuda": tool_version("nvcc", "--version"),
            "mpi": tool_version("mpirun", "--version"),
        },
        "hardware": {
            "cpu": platform.processor() or "unknown",
            "gpu": tool_version(
                "nvidia-smi", "--query-gpu=name", "--format=csv,noheader"
            ),
        },
        "case": None,
    }
    if args.case:
        manifest["case"] = {
            "path": str(args.case),
            "sha256": file_sha256(args.case),
            "size_bytes": args.case.stat().st_size,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
