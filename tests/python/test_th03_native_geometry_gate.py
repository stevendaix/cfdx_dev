"""Regression tests for the fail-closed TH-03 native geometry contract."""

import importlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
native_case_path = importlib.import_module(
    "scripts.require_th03_native_geometry"
).native_case_path
_boundary_patch_integrity_hash = importlib.import_module(
    "cfdx.io.hdf5_writer"
)._boundary_patch_integrity_hash
_boundary_patch_metadata = importlib.import_module(
    "cfdx.io.hdf5_writer"
)._boundary_patch_metadata
make_structured_mesh = importlib.import_module(
    "scripts.generate_th03_openfoam_blockmesh"
).make_structured_mesh


class TestTh03NativeGeometryGate(unittest.TestCase):
    def test_reconstructed_mesh_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "requires native OpenFOAM polyMesh"):
            native_case_path("OpenFOAM v2606 six-block geometry, graded reconstruction")

    def test_missing_provenance_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "requires native OpenFOAM polyMesh"):
            native_case_path("")

    def test_native_prefix_requires_existing_case(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(FileNotFoundError):
            native_case_path(f"OpenFOAM polyMesh, {Path(tmp) / 'missing'}")


if __name__ == "__main__":
    unittest.main()


def test_boundary_patch_integrity_hash_covers_names_types_and_face_ownership():
    mesh = make_structured_mesh("smoke")
    metadata, face_ids, offsets = _boundary_patch_metadata(mesh)
    original = _boundary_patch_integrity_hash(metadata, face_ids, offsets)
    changed_faces = face_ids.copy()
    changed_faces[0], changed_faces[-1] = changed_faces[-1], changed_faces[0]
    assert _boundary_patch_integrity_hash(metadata, changed_faces, offsets) != original
    entries = metadata.split(";")
    first = entries[0].split(":")
    first[3] = str((int(first[3]) + 1) % 8)
    entries[0] = ":".join(first)
    changed_metadata = ";".join(entries)
    assert _boundary_patch_integrity_hash(changed_metadata, face_ids, offsets) != original
