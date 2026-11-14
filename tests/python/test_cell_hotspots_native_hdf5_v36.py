"""V36: compare native HDF5 production OFF/ON residual tracing without meshio.

The quick BFS mesh has a trivial uniform initial state: this certifies only
software paths, transactional trace I/O and diagnostic non-intrusiveness.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import tempfile

from cfdx.cell_convergence import read_cell_trace


def run(solver: str, mesh: str, directory: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [solver, "--mesh", mesh, "--output-dir", str(directory), "--iterations", "2", *extra],
        check=False, capture_output=True, text=True, timeout=75,
    )


def main() -> None:
    solver, mesh = os.environ["CFDX_PRODUCTION_SOLVER"], os.environ["CFDX_PRODUCTION_MESH"]
    if not Path(solver).is_file() or not Path(mesh).is_file():
        raise AssertionError("V36 native HDF5 regression requires a built solver and quick test mesh")
    with tempfile.TemporaryDirectory(prefix="cfdx-v36-cells-") as root:
        root = Path(root)
        off, on = root / "off", root / "on"
        baseline = run(solver, mesh, off)
        recorded = run(solver, mesh, on, "--cell-hotspots", "--cell-hotspots-every", "1", "--cell-map-full")
        assert baseline.returncode == 0, baseline.stderr[-4000:]
        assert recorded.returncode == 0, recorded.stderr[-4000:]
        assert not (off / "cell_convergence.jsonl").exists()
        assert not (off / "cell_convergence_map.bin").exists()
        trace = read_cell_trace(on / "cell_convergence.jsonl")
        assert trace.complete and trace.accepted_iterations == 2
        assert trace.n_cells > 0 and trace.n_faces > 0
        assert [frame.iteration for frame in trace.frames] == [1, 2]
        assert not trace.to_dict()["solver_convergence_certificate"]
        binary = (on / "cell_convergence_map.bin").read_bytes()
        assert binary.startswith(b"CFDXMAP1")
        assert binary[-12:-8] == b"CFND"
        assert int.from_bytes(binary[-8:], "little") == 2
        a = json.loads((off / "convergence.json").read_text())
        b = json.loads((on / "convergence.json").read_text())
        # Diagnostics must not change the nonlinear history, gate, or plan.
        assert a == b, "cell trace changed the underlying FV solution or its convergence evidence"
        assert json.loads((off / "execution.json").read_text())["converged"] is True
        assert json.loads((on / "execution.json").read_text())["converged"] is True
        invalid = run(solver, mesh, root / "bad", "--cell-hotspots-every", "0")
        assert invalid.returncode == 2
        assert "--cell-hotspots-every" in invalid.stderr
    print("V36 NATIVE HDF5: opt-in/off physical evidence identical, maps valid, bounds fail closed")


if __name__ == "__main__":
    main()
