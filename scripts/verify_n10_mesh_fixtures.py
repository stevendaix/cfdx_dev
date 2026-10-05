#!/usr/bin/env python3
"""Verify acquired N10 public mesh fixtures before they can become CI evidence.

This is an integrity gate, not a promotion step. A fixture is VERIFIED only
when the acquired bytes match a non-null manifest SHA-256 and the source is
publicly redistributable. Git-LFS pointer files are explicitly rejected as
mesh payloads because their pointer bytes are not the mesh bytes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.acquire_mesh_fixtures import sha256_tree


LFS_PREFIX = b"version https://git-lfs.github.com/spec/v1\n"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_lfs_pointer(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(len(LFS_PREFIX)) == LFS_PREFIX
    except OSError:
        return False


def parse_manifest(path: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    sources: dict[str, dict[str, str]] = {}
    fixtures: list[dict[str, str]] = []
    section = None
    current: dict[str, str] | None = None

    def finish() -> None:
        nonlocal current
        if current is None:
            return
        if section == "sources":
            sources[current["id"]] = current
        elif section == "fixtures":
            fixtures.append(current)
        current = None

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line in {"sources:", "fixtures:"}:
            finish()
            section = line[:-1]
            continue
        if line.startswith("version:") or line.startswith("schema:"):
            continue
        if line.startswith("- id:"):
            finish()
            current = {"id": line.split(":", 1)[1].strip()}
            continue
        if current is not None and ":" in line:
            key, value = (part.strip() for part in line.split(":", 1))
            current[key] = value
    finish()
    return sources, fixtures


def verify_fixture(
    fixture: dict[str, str], source: dict[str, str], output_root: Path
) -> dict[str, object]:
    fixture_id = fixture["id"]
    result: dict[str, object] = {
        "id": fixture_id,
        "source": fixture.get("source"),
        "status": "BLOCKED",
        "reason": "",
    }

    if fixture.get("source") == "ansys-vmfl":
        result["reason"] = "proprietary ANSYS VMFL is external-only"
        return result

    if source.get("redistribution") not in {"allowed", "verified"}:
        result["reason"] = "source redistribution status is not verified"
        return result

    digest = fixture.get("sha256")
    if digest in {None, "", "null"}:
        result["reason"] = "manifest SHA-256 is not verified"
        return result

    payload = output_root / fixture_id
    if not payload.exists():
        result["reason"] = "acquired fixture is missing"
        return result

    if is_lfs_pointer(payload):
        result["reason"] = "acquired file is a Git-LFS pointer, not mesh payload bytes"
        return result

    if payload.is_file():
        actual = sha256(payload)
    elif payload.is_dir() and fixture.get("acquisition") == "pinned_repository_subtree":
        actual = sha256_tree(payload)
    else:
        result["reason"] = "unsupported acquired fixture payload type"
        return result

    result["sha256"] = actual
    result["expected_sha256"] = digest
    if actual.lower() != digest.lower():
        result["reason"] = "SHA-256 mismatch"
        return result

    result["status"] = "VERIFIED"
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest", type=Path, default=Path("tests/fixtures/mesh_sources.yaml")
    )
    parser.add_argument(
        "--output-root", type=Path, default=Path("build/n10-fixtures")
    )
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()

    sources, fixtures = parse_manifest(args.manifest)
    results = []
    for fixture in fixtures:
        source = sources.get(fixture.get("source", ""), {})
        results.append(verify_fixture(fixture, source, args.output_root))

    verified = [item for item in results if item["status"] == "VERIFIED"]
    report = {
        "campaign": "N10 public mesh fixture verification",
        "status": "PASS" if len(verified) == len(results) else "BLOCKED",
        "fixtures": results,
        "promotion_policy": {
            "require_exact_sha256": True,
            "require_verified_redistribution": True,
            "reject_git_lfs_pointer_as_payload": True,
            "proprietary_ansys_external_only": True,
            "no_automatic_manifest_mutation": True,
        },
    }

    report_path = args.report or args.output_root / "n10_fixture_verification.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
