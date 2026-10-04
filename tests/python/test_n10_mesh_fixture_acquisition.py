from pathlib import Path

import pytest

from scripts.acquire_mesh_fixtures import parse_manifest, parse_sources, sha256


def test_manifest_parsers_preserve_sections():
    manifest = Path("tests/fixtures/mesh_sources.yaml")
    fixtures = parse_manifest(manifest)
    sources = parse_sources(manifest)

    assert len(sources) == 4
    assert len(fixtures) == 5
    assert sources["ansys-vmfl"]["redistribution"] == "prohibited"
    assert fixtures[-1]["source"] == "ansys-vmfl"


def test_sha256_is_byte_exact(tmp_path):
    payload = tmp_path / "mesh"
    payload.write_bytes(b"CFDX-N10\n")
    assert sha256(payload) == "d5c8085b9f4c5b9a2f11f70f22b0e8b3b2b7e42b2b6a6a2c3b9f6d9b1e4e4e6a"
