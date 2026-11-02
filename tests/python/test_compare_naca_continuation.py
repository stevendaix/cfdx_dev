import csv
import tempfile
import unittest
from pathlib import Path

from scripts.compare_naca_continuation import FIELDS, analyze


class ComparisonTest(unittest.TestCase):
    def test_missing_data_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.csv"
            path.write_text("mesh,strategy\nN256,adaptive\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                analyze(path)

    def test_valid_row_preserves_unqualified_flag(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDS)
                writer.writeheader()
                writer.writerow(
                    {
                        "mesh": "N256",
                        "strategy": "adaptive",
                        "converged": "false",
                        "wall_seconds": "10",
                        "true_momentum_residual": "0.1",
                        "true_continuity_residual": "0.1",
                        "mass_imbalance": "0.01",
                        "cl": "0.3",
                        "cd": "0.02",
                        "cp_error": "0.2",
                        "cf_error": "0.3",
                    }
                )
            self.assertFalse(analyze(path)["results"][0]["converged"])


if __name__ == "__main__":
    unittest.main()
