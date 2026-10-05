#!/usr/bin/env python3
"""Acquire N10 external mesh fixtures deterministically from pinned sources.

The manifest is the source of truth. This tool never changes the manifest and
never treats a newly computed digest as verified metadata automatically.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path
from urllib.parse import quote


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


def sha256_tree(path: Path) -> str:
    """Hash a directory deterministically from relative paths and file bytes."""
    digest = hashlib.sha256()
    for child in sorted(p for p in path.rglob("*") if p.is_file()):
        relative = child.relative_to(path).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        with child.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str | None, fixture_id: str, *, directory: bool = False) -> str:
    actual = sha256_tree(path) if directory else sha256(path)
    if expected not in {None, "", "null"} and actual.lower() != expected.lower():
        fail(
            f"fixture {fixture_id}: sha256 mismatch; "
            f"expected={expected} actual={actual}"
        )
    return actual


LFS_POINTER_PREFIX = b"version https://git-lfs.github.com/spec/v1\n"


def is_lfs_pointer(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(len(LFS_POINTER_PREFIX)) == LFS_POINTER_PREFIX


def github_media_url(repository: str, ref: str, relative_path: str) -> str:
    owner_repo = repository.removeprefix("https://github.com/").rstrip("/")
    encoded_path = quote(relative_path, safe="/")
    return f"https://media.githubusercontent.com/media/{owner_repo}/{ref}/{encoded_path}"


def download_url(url: str, destination: Path) -> None:
    with urllib.request.urlopen(url, timeout=30) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out)


def acquire_file(source: dict[str, str], relative_path: str, destination: Path) -> None:
    repository = source.get("repository", "")
    ref = source.get("ref", "")
    if not repository or repository == "null":
        fail(f"source {source['id']}: no public repository for {relative_path}")
    if not ref:
        fail(f"source {source['id']}: immutable ref required for {relative_path}")
    if repository.startswith("https://github.com/") and relative_path:
        owner_repo = repository.removeprefix("https://github.com/").rstrip("/")
        url = f"https://raw.githubusercontent.com/{owner_repo}/{ref}/{quote(relative_path, safe='/')}"
    elif repository.startswith("https://gitlab.") and relative_path:
        url = f"{repository.rstrip('/')}/-/raw/{ref}/{quote(relative_path, safe='/')}"
    else:
        fail(f"source {source['id']}: unsupported repository URL {repository}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        download_url(url, destination)
        if is_lfs_pointer(destination):
            if not repository.startswith("https://github.com/"):
                fail(
                    f"fixture {relative_path}: source returned a Git-LFS pointer, "
                    "but no GitHub LFS media endpoint is available"
                )
            destination.unlink()
            media_url = github_media_url(repository, ref, relative_path)
            download_url(media_url, destination)
            if is_lfs_pointer(destination):
                fail(
                    f"fixture {relative_path}: GitHub LFS media endpoint still returned "
                    "a pointer; actual payload bytes were not acquired"
                )
    except Exception as exc:
        destination.unlink(missing_ok=True)
        if str(exc).startswith("fixture "):
            fail(str(exc))
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
        try:
            subprocess.run(["git", "init", "--quiet", str(checkout)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            subprocess.run(["git", "-C", str(checkout), "remote", "add", "origin", repository], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            subprocess.run(["git", "-C", str(checkout), "fetch", "--filter=blob:none", "--depth", "1", "origin", ref], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            subprocess.run(["git", "-C", str(checkout), "sparse-checkout", "set", "--no-cone", relative_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            subprocess.run(["git", "-C", str(checkout), "checkout", "--detach", "FETCH_HEAD"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
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
        is_directory = fixture.get("acquisition") == "pinned_repository_subtree"
        if is_directory:
            acquire_directory(source, relative, destination)
        else:
            acquire_file(source, relative, destination)
        actual = verify_sha256(destination, fixture.get("sha256"), fixture["id"], directory=is_directory)
        kind = "directory" if is_directory else "file"
        print(f"{fixture['id']}: acquired {kind} {destination} sha256={actual} expected={fixture.get('sha256', 'null')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
