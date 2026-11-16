"""Read-only V77 wall evidence must fail closed on missing/mutated values."""

import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "scripts"))
from bfs_wall_family_preflight_v77 import read_yplus_counts


class TestBfsWallV77(unittest.TestCase):
    def test_real_format_yplus_counts(self):
        line = (
            "BFS_WALL_REGIME_V74 iteration=23 faces=34 wall_resolved=0 buffer=1 "
            "log_layer=11 out_of_range=22 min_yplus_spalding=3.6435727 "
            "max_yplus_spalding=732.29 qualified=0\n"
        )
        with tempfile.TemporaryDirectory() as p:
            path = pathlib.Path(p) / "log"
            path.write_text(line)
            r = read_yplus_counts(path)
        self.assertEqual(r["faces"], 34)
        self.assertEqual(r["out_of_range"], 22)
        self.assertAlmostEqual(r["face_fraction_out_of_range"], 22 / 34)

    def test_refuses_missing_walls(self):
        with tempfile.TemporaryDirectory() as p:
            path = pathlib.Path(p) / "log"
            path.write_text("no evidence")
            with self.assertRaises(ValueError):
                read_yplus_counts(path)

    def test_refuses_false_total(self):
        with tempfile.TemporaryDirectory() as p:
            path = pathlib.Path(p) / "log"
            path.write_text(
                "BFS_WALL_REGIME_V74 faces=34 wall_resolved=1 buffer=1 "
                "log_layer=11 out_of_range=22 min_yplus_spalding=1 max_yplus_spalding=750\n"
            )
            with self.assertRaises(ValueError):
                read_yplus_counts(path)


if __name__ == "__main__":
    unittest.main()
