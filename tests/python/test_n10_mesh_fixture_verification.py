from pathlib import Path

from scripts.verify_n10_mesh_fixtures import is_lfs_pointer, sha256


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
