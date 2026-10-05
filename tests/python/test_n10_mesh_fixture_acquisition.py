from pathlib import Path

import pytest

from scripts.acquire_mesh_fixtures import parse_manifest, parse_sources, sha256, sha256_tree, verify_sha256


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
    assert sha256(payload) == "d1a8147b7ff4ff49fac7ac2d87b71fce68f95266ca9ec5b6193db5b57c8b7c93"


def test_sha256_mismatch_is_rejected(tmp_path):
    payload = tmp_path / "mesh"
    payload.write_bytes(b"CFDX-N10\n")
    with pytest.raises(SystemExit, match="sha256 mismatch"):
        verify_sha256(payload, "0" * 64, "fixture")


def test_sha256_tree_is_path_and_content_deterministic(tmp_path):
    root = tmp_path / "tree"
    (root / "b").mkdir(parents=True)
    (root / "a").mkdir()
    (root / "b" / "x").write_bytes(b"x")
    (root / "a" / "y").write_bytes(b"y")

    digest = sha256_tree(root)
    assert digest == "4cd3cbace896f7ac9fbb224f47ec34a4dd42e10a88b375362f297837bfa73583"
