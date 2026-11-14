"""Regression tests for the fail-closed TH-03 native geometry contract."""
import tempfile
import unittest
from pathlib import Path

import h5py

from scripts.require_th03_native_geometry import native_case_path


class TestTh03NativeGeometryGate(unittest.TestCase):
    def test_reconstructed_mesh_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "requires native OpenFOAM polyMesh"):
            native_case_path("OpenFOAM v2606 six-block geometry, graded reconstruction")

    def test_missing_provenance_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "requires native OpenFOAM polyMesh"):
            native_case_path("")

    def test_native_prefix_requires_existing_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                native_case_path(f"OpenFOAM polyMesh, {Path(tmp) / 'missing'}")


if __name__ == "__main__":
    unittest.main()
