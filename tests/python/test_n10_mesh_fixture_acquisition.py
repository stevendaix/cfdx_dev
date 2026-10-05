from pathlib import Path

import pytest

from scripts.acquire_mesh_fixtures import parse_manifest, parse_sources, sha256, sha256_tree, verify_sha256
from scripts.verify_n10_mesh_fixtures import verify_fixture


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


def test_public_fixture_provenance_is_fully_verified():
    manifest = Path("tests/fixtures/mesh_sources.yaml")
    sources = parse_sources(manifest)
    fixtures = {item["id"]: item for item in parse_manifest(manifest)}

    public_ids = {
        "meshio-su2-square",
        "meshio-gmsh-insulated-2-2",
        "meshio-vtk-unstructured",
        "openfoam-airfoil2d",
    }
    expected = {
        "meshio-su2-square": "58c56bef6d32ae93e4c8ab6fbd5ab92676bd353ea0d160a702fd6bb101b4b16c",
        "meshio-gmsh-insulated-2-2": "e948c04847537cffd8b00fbcc4d5add982d18c76cab2200635d43b84ded08523",
        "meshio-vtk-unstructured": "bba5171c092aa0ed9714ab5fc9f3715b430500554e58e2b8dac3dd9443988414",
        "openfoam-airfoil2d": "c7fe0122a4c8e2b7aaaaa34a7a403246f9dccf2e1b0005d8ba25eecf73b2b533",
    }

    assert public_ids <= fixtures.keys()
    for fixture_id in public_ids:
        fixture = fixtures[fixture_id]
        source = sources[fixture["source"]]
        assert fixture["sha256"] == expected[fixture_id]
        assert fixture["status"] == "verified_reference"
        assert source["license_status"] == "verified"
        assert source["redistribution"] == "verified"

    assert fixtures["ansys-vmfl078"]["sha256"] == "null"
    assert sources["ansys-vmfl"]["redistribution"] == "prohibited"


def test_independent_verifier_checks_directory_digest(tmp_path):
    output_root = tmp_path / "output"
    payload = output_root / "fixture"
    payload.mkdir(parents=True)
    (payload / "a").write_text("one", encoding="utf-8")
    (payload / "b").write_text("two", encoding="utf-8")

    fixture = {
        "id": "fixture",
        "source": "public-source",
        "sha256": "28fe29adbcf0f5e657979954392f18631a1e495dd11c77be6e58e6841a2a2c2d",
        "acquisition": "pinned_repository_subtree",
    }
    source = {"id": "public-source", "redistribution": "verified"}

    result = verify_fixture(fixture, source, output_root)

    assert result["status"] == "VERIFIED"
    assert result["sha256"] == fixture["sha256"]
