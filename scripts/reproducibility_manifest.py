#!/usr/bin/env python3
"""Emit a machine-readable CFDX reproducibility manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import shutil
import datetime
from pathlib import Path


def git_value(root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


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
        return subprocess.check_output(
            [executable, *args], text=True, stderr=subprocess.STDOUT
        ).strip().splitlines()[0]
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    manifest = {
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
        "cfdx_version": "0.7",
        "schema_version": 1,
        "git": {
            "commit": git_value(root, "rev-parse", "HEAD"),
            "dirty": bool(git_value(root, "status", "--porcelain")),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "compiler": {
            "cxx": git_value(root, "config", "--get", "CXX") or "unknown",
            "version": tool_version(git_value(root, "config", "--get", "CXX") or "c++", "--version"),
        },
        "toolchain": {
            "cuda": tool_version("nvcc", "--version"),
            "mpi": tool_version("mpirun", "--version"),
        },
        "hardware": {
            "cpu": platform.processor() or "unknown",
            "gpu": tool_version("nvidia-smi", "--query-gpu=name", "--format=csv,noheader"),
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
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
