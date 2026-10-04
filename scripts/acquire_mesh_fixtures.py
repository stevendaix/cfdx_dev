#!/usr/bin/env python3
"""Acquire N10 external mesh fixtures deterministically from pinned sources.

The manifest is the source of truth. This tool never changes the manifest and
never treats a newly computed digest as verified metadata automatically.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_manifest(path: Path) -> list[dict[str, str]]:
    items: list[dict[str, str]] = []
    section = None
    current = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line in {"sources:", "fixtures:"}:
            section = line[:-1]
            continue
        if line.startswith("version:") or line.startswith("schema:"):
            continue
        if line.startswith("- id:"):
            if section == "fixtures":
                if current is not None:
                    items.append(current)
                current = {"id": line.split(":", 1)[1].strip()}
            continue
        if section == "fixtures" and current is not None and ":" in line:
            key, value = (part.strip() for part in line.split(":", 1))
            current[key] = value
    if current is not None:
        items.append(current)
    return items


def parse_sources(path: Path) -> dict[str, dict[str, str]]:
    sources: dict[str, dict[str, str]] = {}
    section = None
    current = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line in {"sources:", "fixtures:"}:
            if current is not None and section == "sources":
                sources[current["id"]] = current
                current = None
            section = line[:-1]
            continue
        if line.startswith("version:") or line.startswith("schema:"):
            continue
        if section == "sources" and line.startswith("- id:"):
            if current is not None:
                sources[current["id"]] = current
            current = {"id": line.split(":", 1)[1].strip()}
        elif section == "sources" and current is not None and ":" in line:
            key, value = (part.strip() for part in line.split(":", 1))
            current[key] = value
    if current is not None and section == "sources":
        sources[current["id"]] = current
    return sources


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def acquire_file(source: dict[str, str], relative_path: str, destination: Path) -> None:
    repository = source.get("repository", "")
    ref = source.get("ref", "")
    if not repository or repository == "null":
        fail(f"source {source['id']}: no public repository for {relative_path}")
    if repository.startswith("https://github.com/") and relative_path:
        owner_repo = repository.removeprefix("https://github.com/").rstrip("/")
        url = f"https://raw.githubusercontent.com/{owner_repo}/{ref}/{relative_path}"
    elif repository.startswith("https://gitlab."):
        url = f"{repository.rstrip('/')}/-/raw/{ref}/{relative_path}"
    else:
        fail(f"source {source['id']}: unsupported repository URL {repository}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=30) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)
    except Exception as exc:
        destination.unlink(missing_ok=True)
        fail(f"download failed for {relative_path}: {exc}")


def acquire_directory(source: dict[str, str], relative_path: str, destination: Path) -> None:
    repository = source.get("repository", "")
    ref = source.get("ref", "")
    if not repository or repository == "null" or not ref:
        fail(f"source {source['id']}: repository/ref required for directory acquisition")
    if destination.exists():
        shutil.rmtree(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cfdx-n10-fixture-") as tmp:
        checkout = Path(tmp) / "repo"
        cmd = [
            "git", "clone", "--filter=blob:none", "--no-checkout",
            "--depth", "1", "--revision", ref, repository, str(checkout)
        ]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            subprocess.run(
                ["git", "-C", str(checkout), "sparse-checkout", "set", relative_path],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            subprocess.run(
                ["git", "-C", str(checkout), "checkout", ref],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            fail(f"git acquisition failed for {relative_path}: {exc}")
        source_dir = checkout / relative_path
        if not source_dir.is_dir():
            fail(f"source {source['id']}: directory not found at {relative_path}")
        shutil.copytree(source_dir, destination)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("tests/fixtures/mesh_sources.yaml"))
    parser.add_argument("--fixture", action="append", dest="fixtures")
    parser.add_argument("--output-root", type=Path, default=Path("build/n10-fixtures"))
    parser.add_argument("--allow-unverified", action="store_true",
                        help="permit acquisition when manifest sha256 is null; never marks it verified")
    args = parser.parse_args()

    sources = parse_sources(args.manifest)
    fixtures = parse_manifest(args.manifest)
    selected = set(args.fixtures or [item["id"] for item in fixtures])

    for fixture in fixtures:
        if fixture["id"] not in selected:
            continue
        if fixture.get("source") == "ansys-vmfl":
            fail(f"fixture {fixture['id']}: proprietary source requires an authorized external artifact")
        if fixture.get("sha256") in {None, "null"} and not args.allow_unverified:
            fail(f"fixture {fixture['id']}: sha256 is null; use --allow-unverified only for acquisition staging")
        source = sources.get(fixture.get("source", ""))
        if source is None:
            fail(f"fixture {fixture['id']}: unknown source {fixture.get('source')}")
        relative = fixture["path"]
        destination = args.output_root / fixture["id"]
        if fixture.get("acquisition") == "pinned_repository_subtree":
            acquire_directory(source, relative, destination)
        else:
            acquire_file(source, relative, destination)
        if destination.is_file():
            actual = sha256(destination)
            print(f"{fixture['id']}: sha256={actual} expected={fixture.get('sha256', 'null')}")
        else:
            print(f"{fixture['id']}: acquired directory {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
