#!/usr/bin/env python3
"""Validate the N10 external mesh-fixture manifest without fetching external data.

The validator is intentionally dependency-free. It checks the manifest's
integrity boundary: immutable source refs, explicit acquisition policy, valid
SHA-256 values when present, and the rule that mandatory fixtures cannot have
an unverified/null hash.
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

MANIFEST = Path("tests/fixtures/mesh_sources.yaml")
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
HEX40_RE = re.compile(r"^[0-9a-fA-F]{40}$")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_manifest(text: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Parse the deliberately flat, schema-v1 manifest conservatively."""
    sources: list[dict[str, str]] = []
    fixtures: list[dict[str, str]] = []
    section = None
    current = None

    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line == "sources:":
            section = "sources"
            continue
        if line == "fixtures:":
            section = "fixtures"
            continue
        if line.startswith("version:") or line.startswith("schema:"):
            continue
        if line.startswith("- id:"):
            if current is not None:
                (sources if section == "sources" else fixtures).append(current)
            current = {"id": line.split(":", 1)[1].strip()}
            if section not in {"sources", "fixtures"}:
                fail(f"{MANIFEST}:{lineno}: entry outside sources/fixtures")
            continue
        if current is None:
            fail(f"{MANIFEST}:{lineno}: unexpected content")
        if ":" not in line:
            fail(f"{MANIFEST}:{lineno}: expected key: value")
        key, value = (part.strip() for part in line.split(":", 1))
        current[key] = value

    if current is not None:
        (sources if section == "sources" else fixtures).append(current)
    return sources, fixtures


def main() -> int:
    if not MANIFEST.is_file():
        fail(f"missing {MANIFEST}")

    text = MANIFEST.read_text(encoding="utf-8")
    sources, fixtures = parse_manifest(text)

    if "version: 1" not in text or "schema: 1" not in text:
        fail("manifest must declare version: 1 and schema: 1")

    source_ids = {item.get("id") for item in sources}
    if len(source_ids) != len(sources):
        fail("duplicate source id")
    fixture_ids = {item.get("id") for item in fixtures}
    if len(fixture_ids) != len(fixtures):
        fail("duplicate fixture id")

    for source in sources:
        sid = source["id"]
        ref = source.get("ref", "")
        if sid != "ansys-vmfl" and not HEX40_RE.fullmatch(ref):
            fail(f"source {sid}: ref must be an immutable 40-hex commit")
        if sid == "ansys-vmfl" and source.get("redistribution") != "prohibited":
            fail("ANSYS VMFL redistribution must remain prohibited")

    for fixture in fixtures:
        fid = fixture["id"]
        source = fixture.get("source")
        if source not in source_ids:
            fail(f"fixture {fid}: unknown source {source!r}")
        digest = fixture.get("sha256")
        mandatory = fixture.get("mandatory_ci") == "true"
        reference_only = fixture.get("reference_only") == "true"

        if digest not in {None, "null"} and not SHA256_RE.fullmatch(digest):
            fail(f"fixture {fid}: sha256 must be 64 hex characters or null")
        if mandatory and digest in {None, "null"}:
            fail(f"fixture {fid}: mandatory_ci=true requires a verified sha256")
        if mandatory and reference_only:
            fail(f"fixture {fid}: mandatory_ci=true conflicts with reference_only=true")
        if source == "ansys-vmfl" and not reference_only:
            fail(f"fixture {fid}: proprietary ANSYS material must remain reference_only")

    print(f"OK: {MANIFEST} validated ({len(sources)} sources, {len(fixtures)} fixtures)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
