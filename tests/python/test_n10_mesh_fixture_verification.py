from pathlib import Path

from scripts.acquire_mesh_fixtures import github_media_url, sha256_tree
from scripts.verify_n10_mesh_fixtures import is_lfs_pointer, sha256, verify_fixture


def test_lfs_pointer_is_not_payload(tmp_path):
    payload = tmp_path / "mesh"
    payload.write_bytes(
        b"version https://git-lfs.github.com/spec/v1\n"
        b"oid sha256:0123456789abcdef\nsize 12\n"
    )
    assert is_lfs_pointer(payload)


def test_real_payload_is_not_lfs_pointer(tmp_path):
    payload = tmp_path / "mesh"
    payload.write_bytes(b"CFDX-N10\n")
    assert not is_lfs_pointer(payload)
    assert sha256(payload) == (
        "d1a8147b7ff4ff49fac7ac2d87b71fce68f95266ca9ec5b6193db5b57c8b7c93"
    )


def test_github_lfs_media_url_is_pinned_to_ref():
    assert github_media_url(
        "https://github.com/nschloe/meshio",
        "b2ee99842e119901349fdeee06b5bf61e01f450a",
        "tests/meshes/su2/square.su2",
    ) == (
        "https://media.githubusercontent.com/media/nschloe/meshio/"
        "b2ee99842e119901349fdeee06b5bf61e01f450a/tests/meshes/su2/square.su2"
    )


def test_tree_digest_changes_when_payload_changes(tmp_path):
    fixture = tmp_path / "fixture"
    fixture.mkdir()
    (fixture / "a").write_text("one", encoding="utf-8")
    (fixture / "b").write_text("two", encoding="utf-8")
    first = sha256_tree(fixture)

    (fixture / "b").write_text("three", encoding="utf-8")
    second = sha256_tree(fixture)

    assert first == "28fe29adbcf0f5e657979954392f18631a1e495dd11c77be6e58e6841a2a2c2d"
    assert second == "584cb1b22d965daff01fc35dc24c18ddc972a56da7dc925ecbe1258063b3aa00"
    assert first != second


def test_directory_fixture_uses_tree_digest(tmp_path):
    output_root = tmp_path / "output"
    payload = output_root / "openfoam-airfoil2d"
    payload.mkdir(parents=True)
    (payload / "points").write_text("points", encoding="utf-8")
    (payload / "faces").write_text("faces", encoding="utf-8")
    expected = sha256_tree(payload)

    fixture = {
        "id": "openfoam-airfoil2d",
        "source": "openfoam-dev",
        "sha256": expected,
        "acquisition": "pinned_repository_subtree",
    }
    source = {"id": "openfoam-dev", "redistribution": "verified"}

    result = verify_fixture(fixture, source, output_root)

    assert result["status"] == "VERIFIED"
    assert result["sha256"] == expected
    assert result["expected_sha256"] == expected
