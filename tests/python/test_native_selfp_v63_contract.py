import importlib.util
import tempfile
import unittest
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[2] / "scripts/qualify_native_selfp_v63.py"
spec = importlib.util.spec_from_file_location("qualify_native_selfp_v63", SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class V63Contract(unittest.TestCase):
    def test_six_methods_are_distinct_and_complete(self):
        self.assertEqual(len(module.METHODS), 6)
        self.assertEqual(len(set(module.METHODS)), 6)
        self.assertIn("pcd", module.METHODS)
        self.assertNotIn("selfp", module.METHODS)

    def test_missing_captures_rejected(self):
        with tempfile.TemporaryDirectory() as path:
            binary = Path(path) / "fake_replay"
            binary.write_text("#!/bin/sh\nexit 0\n")
            with self.assertRaises(ValueError):
                module.audit(binary, Path(path))

    def test_identical_operators_required(self):
        with tempfile.TemporaryDirectory() as path:
            directory = Path(path)
            binary = directory / "fake_replay"
            binary.write_text("#!/bin/sh\nexit 0\n")
            for method in module.METHODS:
                for kind in ("A", "b", "x0"):
                    (directory / f"{method}_{kind}.mtx").write_text("data")
            (directory / "pcd_A.mtx").write_text("different")
            with self.assertRaisesRegex(ValueError, "input mismatch"):
                module.audit(binary, directory)

    def test_false_success_rejected(self):
        with tempfile.TemporaryDirectory() as path:
            directory = Path(path)
            binary = directory / "fake_replay"
            binary.write_text(
                '#!/bin/sh\necho "V63_SELFP_REPLAY n=608 nnz=7080 schur_nnz=1212 iterations=3 true_relative=1e-3 claimed_relative=1e-12 status=0 physical_qualified=0"\n'
            )
            binary.chmod(0o755)
            for method in module.METHODS:
                for kind in ("A", "b", "x0"):
                    (directory / f"{method}_{kind}.mtx").write_text("data")
            with self.assertRaisesRegex(ValueError, "not certified"):
                module.audit(binary, directory)


if __name__ == "__main__":
    unittest.main()
