import tempfile
import unittest
from pathlib import Path

from scripts.plot_bfs_pressure_v79 import extract_native_oracle, run


def example(csr, face, commit="0"):
    return (
        "PRESSURE_SPLIT_ORACLE_V79 iter=1 corrector=0 "
        f"csr_relative={csr} double_only_relative=3e-10 "
        f"physical_face_relative={face} csr_closed={int(csr <= 1e-11)} "
        f"physical_face_closed={int(face <= 1e-11)} requested_tolerance=1e-11 "
        f"production_commit_authorized={commit} read_only=1\n"
    )


class PlotV79(unittest.TestCase):
    def test_valid_readonly_comparison(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)
            (p / "a.log").write_text(example(1e-13, 3e-10))
            (p / "b.log").write_text(example(3e-10, 8e-14))
            r = run(p / "a.log", p / "b.log", p / "out")
            self.assertFalse(r["physical_validation_authorized"])
            self.assertGreater(r["physical_face_improvement_factor"], 1000)
            self.assertTrue((p / "out" / "01_csr_vs_rhie_chow_v79.png").is_file())

    def test_reject_fictitious_commit(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "bad.log"
            p.write_text(example(1e-13, 3e-10, "1"))
            with self.assertRaises(ValueError):
                extract_native_oracle(p)

    def test_reject_incoherent_gates(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "bad.log"
            p.write_text(example(1e-13, 3e-10).replace("csr_closed=1", "csr_closed=0"))
            with self.assertRaises(ValueError):
                extract_native_oracle(p)


if __name__ == "__main__":
    unittest.main()
