"""V36: compare native HDF5 production OFF/ON residual tracing without meshio.

The quick BFS mesh has a trivial uniform initial state: this certifies only
software paths, transactional trace I/O and diagnostic non-intrusiveness.
"""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path

from cfdx.cell_convergence import read_cell_trace


def run(
    solver: str, mesh: str, directory: Path, *extra: str
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            solver,
            "--mesh",
            mesh,
            "--output-dir",
            str(directory),
            "--iterations",
            "2",
            *extra,
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=75,
    )


def main() -> None:
    solver, mesh = (
        os.environ["CFDX_PRODUCTION_SOLVER"],
        os.environ["CFDX_PRODUCTION_MESH"],
    )
    if not Path(solver).is_file() or not Path(mesh).is_file():
        raise AssertionError(
            "V36 native HDF5 regression requires a built solver and quick test mesh"
        )
    with tempfile.TemporaryDirectory(prefix="cfdx-v36-cells-") as root:
        root = Path(root)
        off, on, budget = root / "off", root / "on", root / "budget"
        budget_top = root / "budget_top"
        baseline = run(solver, mesh, off)
        recorded = run(
            solver,
            mesh,
            on,
            "--cell-hotspots",
            "--cell-hotspots-every",
            "1",
            "--cell-map-full",
        )
        budgeted = run(
            solver,
            mesh,
            budget,
            "--momentum-budget-csv",
            "momentum_face_budget.csv",
            "--momentum-budget-every",
            "1",
            "--momentum-budget-all-cells",
        )
        budgeted_top = run(
            solver,
            mesh,
            budget_top,
            "--momentum-budget-csv",
            "momentum_face_budget.csv",
            "--momentum-budget-every",
            "1",
        )
        assert baseline.returncode == 0, baseline.stderr[-4000:]
        assert recorded.returncode == 0, recorded.stderr[-4000:]
        assert budgeted.returncode == 0, budgeted.stderr[-4000:]
        assert budgeted_top.returncode == 0, budgeted_top.stderr[-4000:]
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
        c = json.loads((budget / "convergence.json").read_text())
        d = json.loads((budget_top / "convergence.json").read_text())
        # Neither diagnostic may change the nonlinear history, gate, or plan.
        assert a == b == c == d, (
            "diagnostic output changed the underlying FV solution or convergence evidence"
        )
        assert json.loads((off / "execution.json").read_text())["converged"] is True
        assert json.loads((on / "execution.json").read_text())["converged"] is True
        assert json.loads((budget / "execution.json").read_text())["converged"] is True

        budget_path = budget / "momentum_face_budget.csv"
        with budget_path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert rows, "accepted momentum budget callback must emit CSV rows"
        expected_columns = {
            "iteration",
            "component",
            "rank",
            "selection_reason",
            "cell",
            "x",
            "y",
            "z",
            "convection",
            "viscous",
            "transpose_viscous",
            "pressure",
            "body_force",
            "bounded_divergence",
            "mass_divergence",
            "collocated_difference",
            "physical_residual",
            "bounded_operator_residual",
            "worst_absolute_residual",
            "independent_face_residual_linf",
            "independent_face_residual_difference",
            "temporal",
            "algebraic_residual",
            "matrix_face_difference",
            "matrix_face_scale",
            "matrix_closure_linf",
            "matrix_closure_relative_linf",
            "matrix_closure_checked",
            "matrix_closure_consistent",
            "matrix_closure_worst_cell",
        }
        assert set(rows[0]) == expected_columns
        groups: dict[tuple[int, int], set[int]] = {}
        for row in rows:
            iteration = int(row["iteration"])
            component = int(row["component"])
            cell = int(row["cell"])
            assert iteration in (1, 2)
            assert component in (0, 1, 2)
            assert row["selection_reason"] == "all_cells"
            assert 0 <= cell < trace.n_cells
            groups.setdefault((iteration, component), set()).add(cell)
            scale = float(row["matrix_face_scale"])
            difference = float(row["matrix_face_difference"])
            independent_face_linf = float(row["independent_face_residual_linf"])
            independent_face_difference = float(
                row["independent_face_residual_difference"]
            )
            worst_absolute = float(row["worst_absolute_residual"])
            assert math.isfinite(scale) and math.isfinite(difference)
            assert math.isfinite(independent_face_linf)
            assert math.isfinite(independent_face_difference)
            assert math.isclose(
                independent_face_difference,
                abs(worst_absolute - independent_face_linf),
                rel_tol=1e-14,
                abs_tol=1e-30,
            )
            assert independent_face_difference <= 1e-10 * max(
                1.0, abs(worst_absolute), abs(independent_face_linf)
            )
            assert row["matrix_closure_checked"] == "1"
            assert row["matrix_closure_consistent"] == "1"
            assert abs(difference) <= 2e-10 * max(scale, 1e-15)
        assert set(groups) == {
            (iteration, component) for iteration in (1, 2) for component in (0, 1, 2)
        }
        assert all(cells == set(range(trace.n_cells)) for cells in groups.values())

        # Compact mode must retain the largest physical residual cells AND
        # the independently selected worst face/CSR closure cell.
        top_path = budget_top / "momentum_face_budget.csv"
        with top_path.open(newline="") as stream:
            top_rows = list(csv.DictReader(stream))
        top_groups: dict[tuple[int, int], list[dict[str, str]]] = {}
        for row in top_rows:
            key = (int(row["iteration"]), int(row["component"]))
            top_groups.setdefault(key, []).append(row)
        assert set(top_groups) == set(groups)
        for selected_rows in top_groups.values():
            assert 1 <= len(selected_rows) <= 13
            worst_closure_cell = int(selected_rows[0]["matrix_closure_worst_cell"])
            selected_cells = {int(row["cell"]) for row in selected_rows}
            assert worst_closure_cell in selected_cells
            assert all(
                row["selection_reason"]
                in {"physical_residual_top12", "matrix_closure_worst"}
                for row in selected_rows
            )
            assert all(
                row["cell"] == row["matrix_closure_worst_cell"]
                for row in selected_rows
                if row["selection_reason"] == "matrix_closure_worst"
            )

        invalid = run(solver, mesh, root / "bad", "--cell-hotspots-every", "0")
        assert invalid.returncode == 2
        assert "--cell-hotspots-every" in invalid.stderr
        invalid_budget = run(
            solver,
            mesh,
            root / "bad_budget",
            "--momentum-budget-all-cells",
        )
        assert invalid_budget.returncode == 2
        assert "--momentum-budget-csv" in invalid_budget.stderr
    print(
        "V36 NATIVE HDF5: diagnostics are nonintrusive; "
        "cell maps and momentum face/CSR budget valid"
    )


if __name__ == "__main__":
    main()
