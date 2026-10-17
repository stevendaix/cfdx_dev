from __future__ import annotations

import json
import os
import math
import re
import subprocess
from pathlib import Path


def test_n10_verified_fixture_reaches_production_solver(tmp_path: Path) -> None:
    """Exercise a verified public fixture through the production solver.

    The mesh-to-CFDX conversion is deliberately performed by the production
    C++ mesh importer, not by Python format auto-detection. This keeps the N10
    fixture path aligned with the independently qualified importer and avoids
    introducing generic VTK support into the Python converter as part of this
    increment.

    This is numerical-support smoke evidence, not physical validation.
    """
    fixture_root_value = os.environ.get("CFDX_N10_FIXTURE_ROOT")
    solver_value = os.environ.get("CFDX_PRODUCTION_SOLVER")
    case_builder_value = os.environ.get("CFDX_N10_PRODUCTION_CASE_BUILDER")
    evidence_value = os.environ.get("CFDX_N10_SOLVER_EVIDENCE")
    if not fixture_root_value or not solver_value or not case_builder_value:
        import pytest

        pytest.skip(
            "N10 production-solver fixture environment is provided by the dedicated workflow"
        )

    fixture_root = Path(fixture_root_value)
    solver = Path(solver_value)
    case_builder = Path(case_builder_value)
    source = fixture_root / "meshio-vtk-unstructured"

    assert source.is_file()
    assert solver.is_file()
    assert case_builder.is_file()

    vtk_source = tmp_path / "06_unstructured.vtk"
    vtk_source.write_bytes(source.read_bytes())

    case_path = tmp_path / "unstructured.cfdx.h5"
    build_case = subprocess.run(
        [str(case_builder), str(vtk_source), str(case_path)],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
        env=os.environ.copy(),
    )
    build_diagnostics = build_case.stdout + build_case.stderr
    assert build_case.returncode == 0, build_diagnostics[-12000:]
    assert case_path.is_file()

    output_dir = tmp_path / "run"
    result = subprocess.run(
        [
            str(solver),
            "--mesh",
            str(case_path),
            "--output-dir",
            str(output_dir),
            "--iterations",
            "20",
        ],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
        env=os.environ.copy(),
    )
    diagnostics = result.stdout + result.stderr
    restart_path = output_dir / "restart.dat"
    assert result.returncode == 0, diagnostics[-12000:]

    # Parse named fields rather than relying on their output order. The solver's
    # trace was extended by N10 to expose separate RHS- and first-iteration-
    # reference-normalized momentum residuals; accept the old field name too
    # for compatibility with traces produced before that change.
    scalar_pattern = re.compile(
        r"(?<!\S)([A-Za-z_]\w*)=([-+0-9.eE]+)(?=\s|$)"
    )
    history: list[dict[str, float | int | list[float]]] = []
    for line in diagnostics.splitlines():
        stripped = line.strip()
        if not stripped.startswith("INCOMPRESSIBLE_ITER "):
            continue
        fields = dict(scalar_pattern.findall(stripped))
        momentum_key = (
            "momentum_residual_reference_relative"
            if "momentum_residual_reference_relative" in fields
            else "momentum_residual_relative"
        )
        required = {
            "iter",
            momentum_key,
            "continuity_normalized",
            "velocity_change_inf",
            "pressure_change_inf",
            "alpha_u",
            "alpha_p",
            "pressure_linear_iterations",
            "pressure_relative_residual",
            "momentum_linear_iterations",
            "flux_velocity_mismatch_linf",
        }
        missing = sorted(required.difference(fields))
        assert not missing, (
            f"Malformed INCOMPRESSIBLE_ITER trace; missing fields {missing}: {stripped}"
        )

        def parse_tuple(name: str) -> list[float]:
            tuple_pattern = re.compile(
                rf"(?<!\S){re.escape(name)}=\(\s*([-+0-9.eE]+)\s*,\s*"
                rf"([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)"
            )
            match = tuple_pattern.search(stripped)
            assert match is not None, (
                f"Malformed INCOMPRESSIBLE_ITER tuple {name}: {stripped}"
            )
            return [float(match.group(index)) for index in range(1, 4)]

        row: dict[str, float | int | list[float]] = {
            "iteration": int(fields["iter"]),
            "momentum_residual_relative": float(fields[momentum_key]),
            "momentum_residual_rhs_relative": float(
                fields.get("momentum_residual_rhs_relative", fields[momentum_key])
            ),
            "continuity_normalized": float(fields["continuity_normalized"]),
            "velocity_change_inf": float(fields["velocity_change_inf"]),
            "pressure_change_inf": float(fields["pressure_change_inf"]),
            "alpha_u": float(fields["alpha_u"]),
            "alpha_p": float(fields["alpha_p"]),
            "pressure_linear_iterations": int(fields["pressure_linear_iterations"]),
            "pressure_relative_residual": float(fields["pressure_relative_residual"]),
            "momentum_linear_iterations": int(fields["momentum_linear_iterations"]),
            "flux_velocity_mismatch_linf": float(fields["flux_velocity_mismatch_linf"]),
            "momentum_rhs_scale_components": parse_tuple("momentum_rhs_scale_components"),
            "momentum_rhs_relative_components": parse_tuple(
                "momentum_rhs_relative_components"
            ),
        }
        numeric_values = [
            float(value)
            for value in row.values()
            if not isinstance(value, list)
            and not isinstance(value, int)
        ]
        numeric_values.extend(
            component
            for value in row.values()
            if isinstance(value, list)
            for component in value
        )
        assert all(math.isfinite(value) for value in numeric_values)
        history.append(row)

    assert history, diagnostics[-12000:]
    assert [int(row["iteration"]) for row in history] == sorted(
        int(row["iteration"]) for row in history
    )
    assert "Resolved numerical selections:" in diagnostics
    assert "Converged YES" in diagnostics
    assert restart_path.is_file()

    evidence = {
        "fixture": "meshio-vtk-unstructured",
        "conversion_path": "production_cpp_mesh_importer",
        "case": str(case_path),
        "production_solver": str(solver),
        "execution_status": "PASS",
        "numerical_support": "execution_and_convergence_smoke_only",
        "physical_validation": "NOT_CLAIMED",
        "restart_artifact": True,
        "iteration_budget": 20,
        "quantitative_solver_evidence": {
            "trace_count": len(history),
            "iterations": history,
            "final_iteration": int(history[-1]["iteration"]),
            "max_momentum_residual_relative": max(
                float(row["momentum_residual_relative"]) for row in history
            ),
            "max_continuity_normalized": max(
                float(row["continuity_normalized"]) for row in history
            ),
            "final_momentum_residual_relative": float(
                history[-1]["momentum_residual_relative"]
            ),
            "final_continuity_normalized": float(history[-1]["continuity_normalized"]),
            "final_velocity_change_inf": float(history[-1]["velocity_change_inf"]),
            "final_pressure_change_inf": float(history[-1]["pressure_change_inf"]),
            "final_flux_velocity_mismatch_linf": float(
                history[-1]["flux_velocity_mismatch_linf"]
            ),
        },
        "case_builder_output_tail": build_diagnostics[-4000:],
        "solver_output_tail": diagnostics[-4000:],
    }
    report = (
        Path(evidence_value)
        if evidence_value
        else output_dir / "n10_numerical_support_smoke.json"
    )
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        test_n10_verified_fixture_reaches_production_solver(Path(directory))
