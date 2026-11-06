"""Self-tests for the N15 benchmark contracts."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts" / "n15_benchmark.py"


def run(*args):
    return subprocess.run(
        [sys.executable, str(RUNNER), *args],
        text=True,
        capture_output=True,
        check=False,
    )


class N15Tests(unittest.TestCase):
    def test_run_and_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            case = d / "case.json"
            result = d / "result.json"
            case.write_text(
                json.dumps(
                    {
                        "id": "n15-selftest",
                        "mesh": {"cells": 16},
                        "numerics": {"solver": "test"},
                    }
                )
            )
            p = run(
                "run",
                "--case",
                str(case),
                "--repo",
                str(ROOT),
                "--output",
                str(result),
                "--",
                sys.executable,
                "-c",
                "print('ok')",
            )
            self.assertEqual(p.returncode, 0, p.stderr)
            data = json.loads(result.read_text())
            self.assertEqual(data["schema_version"], "1.0")
            self.assertEqual(data["status"], "PASS")
            self.assertTrue(data["provenance"]["git_sha"])
            self.assertTrue(data["hardware"]["logical_cpus"])
            if sys.platform.startswith("linux") or sys.platform == "darwin":
                self.assertGreater(data["performance"]["peak_rss_bytes"], 0)
                self.assertEqual(data["performance"]["peak_rss_source"], "wait4-child")

    def test_scale_reports_peak_rss_per_child(self):
        if not hasattr(os, "wait4") or not (
            sys.platform.startswith("linux") or sys.platform == "darwin"
        ):
            self.skipTest("per-child RSS requires supported POSIX wait4 semantics")
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            case = d / "case.json"
            result = d / "scale.json"
            case.write_text(
                json.dumps({"id": "n15-rss-isolation", "mesh": {"cells_per_rank": 16}})
            )
            child = "import os,time; data=bytearray(80*1024*1024) if os.environ['CFDX_N15_MPI_RANKS']=='1' else bytearray(1024); time.sleep(.15)"
            p = run(
                "scale",
                "--case",
                str(case),
                "--repo",
                str(ROOT),
                "--kind",
                "strong",
                "--values",
                "1,2",
                "--output",
                str(result),
                "--",
                sys.executable,
                "-c",
                child,
            )
            self.assertEqual(p.returncode, 0, p.stderr)
            runs = json.loads(result.read_text())["runs"]
            first = runs[0]["performance"]["peak_rss_bytes"]
            second = runs[1]["performance"]["peak_rss_bytes"]
            self.assertGreater(first, second * 2)
            self.assertEqual(runs[1]["performance"]["peak_rss_source"], "wait4-child")

    def test_regression(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            b = d / "b.json"
            c = d / "c.json"
            o = d / "o.json"
            b.write_text(
                json.dumps(
                    {"performance": {"total": 10.0}, "provenance": {"git_sha": "base"}}
                )
            )
            c.write_text(
                json.dumps(
                    {"performance": {"total": 12.0}, "provenance": {"git_sha": "head"}}
                )
            )
            p = run(
                "regression",
                "--baseline",
                str(b),
                "--candidate",
                str(c),
                "--threshold-percent",
                "10",
                "--output",
                str(o),
            )
            self.assertEqual(p.returncode, 1)
            r = json.loads(o.read_text())
            self.assertTrue(r["regression"])
            self.assertEqual(r["relative_change_percent"], 20.0)


if __name__ == "__main__":
    unittest.main()
