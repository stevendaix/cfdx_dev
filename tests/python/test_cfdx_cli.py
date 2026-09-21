#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
CLI = ROOT / "scripts" / "cfdx.py"


def make_case(path: Path) -> None:
    with h5py.File(path, "w") as h:
        h.attrs["schema_version"] = 1
        h.create_dataset("points", data=np.zeros((0, 3)))
        h.create_dataset("face_vertices", data=np.zeros(0, dtype=np.uint64))
        h.create_dataset("face_offsets", data=np.zeros(1, dtype=np.uint64))
        h.create_dataset("owner", data=np.zeros(0, dtype=np.uint64))
        h.create_dataset("neighbour", data=np.zeros(0, dtype=np.int64))
        h.create_dataset("cell_faces", data=np.zeros(0, dtype=np.uint64))
        h.create_dataset("cell_offsets", data=np.zeros(1, dtype=np.uint64))


def main() -> int:
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        case = Path(td) / "case.cfdx.h5"
        make_case(case)
        check = subprocess.run(
            [sys.executable, str(CLI), "check", str(case)],
            capture_output=True, text=True, check=False,
        )
        if check.returncode != 0:
            raise SystemExit(check.stderr)
        info = subprocess.run(
            [sys.executable, str(CLI), "info", str(case)],
            capture_output=True, text=True, check=False,
        )
        if info.returncode != 0 or '"dataset_count": 7' not in info.stdout:
            raise SystemExit(info.stdout)
    print("CFDX_CLI_TEST: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
