"""GCI must never turn partial or nonmonotonic CFD data into a certificate."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path


import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from naca_surface_force_oracle import COLUMNS
    sys.path.append(str(SCRIPTS))

from naca_surface_force_oracle import COLUMNS


def load(name: str):
    source = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GCI = load("naca_dual_family_gci_report")
RUNNER = load("run_naca_gci_case")
QUALIFY = load("qualify_naca_grid")

analyze = GCI.analyze
check_wall_shape = GCI.check_wall_shape
assumed_order_two_grid = GCI.assumed_order_two_grid
qualify = QUALIFY.qualify
run_case = RUNNER.run_case


def quality_rows():
    return [
        {
            "n_surface": n,
            "n_radial": 24 * n // 64,
            "radius": 20.0,
            "radial_power": 2.0,
            "radial_refinement_strategy": "nested_base_24",
            "inverted_or_degenerate_cells": 0,
            "minimum_first_layer_spacing": 0.032 / n * 64,
            "mean_first_layer_spacing": 0.034 / n * 64,
            "maximum_cell_edge_aspect_ratio": 10.0,
            "maximum_first_layer_spacing": 0.037 / n * 64,
            "minimum_absolute_cell_area": 0.001 / n,
            "minimum_radial_growth_ratio": 1.0,
            "maximum_radial_growth_ratio": 1.2,
            "maximum_nonorthogonality_degrees": 40.0,
            "p95_nonorthogonality_degrees": 20.0,
            "wall_edge_aspect_ratio_p95": 2.1,
            "wall_edge_aspect_ratio_p50": 1.8,
            "wall_edge_aspect_ratio_max": 10.0,
        }
        for n in (64, 128, 256)
    ]


def test_gci_separately_recovers_known_power_laws():
    cdp_inf, cdv_inf = 0.041, 0.078
    a_p, a_v = 0.37, 1.3
    p_p, p_v = 1.7, 1.2
    rows = []
    for n in (64, 128, 256):
        h = 1 / n
        cdp = cdp_inf + a_p * h**p_p
        cdv = cdv_inf + a_v * h**p_v
        rows.append(
            {
                "case": f"N{n}",
                "h": str(h),
                "converged": "yes",
                "res_mom": "1e-11",
                "res_mass": "1e-11",
                "cd": str(cdp + cdv),
                "cdp": str(cdp),
                "cdv": str(cdv),
                "cl": "0",
                "cm": "0",
            }
        )
    gci = qualify(rows)
    for key, exact, order in (("cdp", cdp_inf, p_p), ("cdv", cdv_inf, p_v)):
        assert gci["quantities"][key]["status"] == "asymptotic_candidate"
        assert gci["quantities"][key]["extrapolated"] == pytest.approx(exact, abs=1e-10)
        assert gci["quantities"][key]["order"] == pytest.approx(order, abs=1e-8)
    assert gci["quantities"]["cd"]["extrapolated"] == pytest.approx(
        cdp_inf + cdv_inf, abs=0.001
    )


def test_total_cd_monotone_is_not_sufficient_for_viscous_gci():
    rows = []
    for n, cd, cdp in ((64, 0.14, 0.042), (128, 0.13, 0.041), (256, 0.125, 0.0421)):
        rows.append(
            {
                "case": f"N{n}",
                "h": str(1 / n),
                "converged": "yes",
                "res_mom": "1e-11",
                "res_mass": "1e-11",
                "cd": str(cd),
                "cdp": str(cdp),
                "cdv": str(cd - cdp),
                "cl": "0",
                "cm": "0",
            }
        )
    gci = qualify(rows)
    assert gci["quantities"]["cd"]["status"] == "asymptotic_candidate"
    assert gci["quantities"]["cdp"]["status"] == "nonmonotone_or_plateau"


def test_bulk_aspect_similarity_vs_maximum_anisotropy():
    rows = quality_rows()
    assert check_wall_shape(rows)["status"] == "REPEATABLE_BULK_WALL_ASPECT"
    rows[2]["wall_edge_aspect_ratio_p95"] = 7.0
    assert check_wall_shape(rows)["status"] == "VARIABLE_BULK_WALL_ASPECT"


def test_missing_cfd_run_refuses_gci(tmp_path):
    for name in ("quasi_isotropic", "controlled_anisotropic", "anisotropic"):
        p = tmp_path / "families" / name
        p.mkdir(parents=True)
        (p / "quality.json").write_text(json.dumps(quality_rows()))
    result = analyze(tmp_path / "families", tmp_path / "runs", Path(sys.executable))
    assert result["status"] == "BLOCKED"
    assert all(f["status"] == "BLOCKED" for f in result["families"].values())
    assert all("gci" not in f for f in result["families"].values())


def test_two_grid_assumed_gci_requires_verified_physical_cases():
    a = {"N": 64, "converged": True, "cd": 0.11, "cdp": 0.04, "cdv": 0.07}
    b = {"N": 128, "converged": True, "cd": 0.12, "cdp": 0.039, "cdv": 0.081}
    provisional = assumed_order_two_grid(a, b)
    assert provisional["status"] == "ASSUMED_ORDER_ONLY_NOT_QUALIFIED"
    assert provisional["components"]["cdp"][
        "gci_relative_conditional"
    ] == pytest.approx(0.001 / 0.039)
    assert provisional["components"]["cdv"][
        "gci_relative_conditional"
    ] == pytest.approx(0.011 / 0.081)
    assert provisional["components"]["cd"][
        "hypothetical_order_based_limit"
    ] == pytest.approx(0.1233333333333)
    b["converged"] = False
    with pytest.raises(ValueError, match="physically verified"):
        assumed_order_two_grid(a, b)



def write_synthetic_raw_wall_faces(path: Path, cdp: float, cdv: float) -> None:
    """Closed four-sided, 32-face synthetic oracle; NEVER physical NACA data."""
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=COLUMNS)
        writer.writeheader()
        for wall_id, wall in enumerate(("left", "right", "bottom", "top")):
            for index in range(8):
                t = (index + 0.5) / 8.0
                x, y, nx, ny = {
                    "left": (0.0, t - 0.5, -1.0, 0.0),
                    "right": (1.0, t - 0.5, 1.0, 0.0),
                    "bottom": (t, -0.5, 0.0, -1.0),
                    "top": (t, 0.5, 0.0, 1.0),
                }[wall]
                writer.writerow(
                    {
                        "face": wall_id * 8 + index,
                        "x": x,
                        "y": y,
                        "z": 0.5,
                        "nx": nx,
                        "ny": ny,
                        "nz": 0.0,
                        "area": 0.125,
                        "pressure": cdp * x,
                        "tau_n_x": -cdv / 4.0,
                        "tau_n_y": 0.0,
                        "tau_n_z": 0.0,
                        "dynamic_pressure": 1.0,
                        "reference_area": 1.0,
                        "reference_chord": 1.0,
                        "x_ref": 0.25,
                        "y_ref": 0.0,
                    }
                )


@pytest.mark.parametrize(
    "schemes,expected",
    [
        (("TVD_VANLEER", "TVD_VANLEER"), "CROSS_MESH_GCI_AGREEMENT_CANDIDATE"),
        (
            ("MUSCL_CELL_LIMITED", "MUSCL_CELL_LIMITED"),
            "CROSS_MESH_GCI_AGREEMENT_CANDIDATE",
        ),
        (("TVD_VANLEER", "MUSCL_CELL_LIMITED"), "BLOCKED"),
    ],
)
def test_full_synthetic_three_grid_chain_exercises_cross_family_gci(
    tmp_path, schemes, expected
):
def write_synthetic_raw_wall_faces(path: Path, cdp: float, cdv: float) -> None:
    """Closed four-sided, 32-face synthetic oracle; NEVER physical NACA data."""
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=COLUMNS)
        writer.writeheader()
        for wall_id, wall in enumerate(("left", "right", "bottom", "top")):
            for index in range(8):
                t = (index + 0.5) / 8.0
                x, y, nx, ny = {
                    "left": (0.0, t - 0.5, -1.0, 0.0),
                    "right": (1.0, t - 0.5, 1.0, 0.0),
                    "bottom": (t, -0.5, 0.0, -1.0),
                    "top": (t, 0.5, 0.0, 1.0),
                }[wall]
                writer.writerow(
                    {
                        "face": wall_id * 8 + index,
                        "x": x,
                        "y": y,
                        "z": 0.5,
                        "nx": nx,
                        "ny": ny,
                        "nz": 0.0,
                        "area": 0.125,
                        "pressure": cdp * x,
                        "tau_n_x": -cdv / 4.0,
                        "tau_n_y": 0.0,
                        "tau_n_z": 0.0,
                        "dynamic_pressure": 1.0,
                        "reference_area": 1.0,
                        "reference_chord": 1.0,
                        "x_ref": 0.25,
                        "y_ref": 0.0,
                    }
                )


@pytest.mark.parametrize(
    "schemes,expected",
    [
        (("TVD_VANLEER", "TVD_VANLEER"), "CROSS_MESH_GCI_AGREEMENT_CANDIDATE"),
        (
            ("MUSCL_CELL_LIMITED", "MUSCL_CELL_LIMITED"),
            "CROSS_MESH_GCI_AGREEMENT_CANDIDATE",
        ),
        (("TVD_VANLEER", "MUSCL_CELL_LIMITED"), "BLOCKED"),
    ],
)
def test_full_synthetic_three_grid_chain_exercises_cross_family_gci(
    tmp_path, schemes, expected
):
    """Synthetic algebraic values exercise gates; not physical NACA evidence."""
    names = ("quasi_isotropic", "anisotropic")
    for name, scheme in zip(names, schemes):
    for name, scheme in zip(names, schemes):
        folder = tmp_path / "families" / name
        folder.mkdir(parents=True)
        rows = quality_rows()
        for row in rows:
            row["surface_parameterization"] = (
                "arc_length" if name == names[0] else "cosine"
            )
        (folder / "quality.json").write_text(json.dumps(rows))
        for n in (64, 128, 256):
            run = tmp_path / "runs" / name / f"n{n}"
            run.mkdir(parents=True)
            mesh = folder / f"naca0012_{name}_n{n}.h5"
            h = 1.0 / n
            offset = 0.0 if name == names[0] else 2e-7
            pressure = 0.04 + offset + 0.3 * h**1.7
            viscous = 0.08 + offset + 1.0 * h**1.4
            (run / "run_summary.json").write_text(
                json.dumps(
                    {
                        "mesh": str(mesh.resolve()),
                        "scheme": scheme,
                        "algorithm": "SIMPLE",
                        "linear_tolerance": 1e-12,
                        "wall_traction": "FACE_WLS",
                        "nonlinear_tolerance": 1e-8,
                        "exit_code": 0,
                        "pressure_boundary_policy": "ZERO_NORMAL_GRADIENT",
                        "kinematic_face_flux_certified": True,
                    }
                )
            )
            with (run / "force.csv").open("w", newline="") as file:
                writer = csv.DictWriter(
                    file,
                    fieldnames=(
                        "case",
                        "h",
                        "converged",
                        "res_mom",
                        "res_mass",
                        "cd",
                        "cdp",
                        "cdv",
                        "cl",
                        "cm",
                    ),
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "case": f"N{n}",
                        "h": h,
                        "converged": "yes",
                        "res_mom": "1e-11",
                        "res_mass": "2e-11",
                        "cd": pressure + viscous,
                        "cdp": pressure,
                        "cdv": viscous,
                        "cl": 0,
                        "cm": 0,
                    }
                )
            write_synthetic_raw_wall_faces(
                run / "raw_wall_faces.csv", pressure, viscous
            )
            (run / "history.csv").write_text(
                "iteration,momentum_rhs_relative,continuity_normalized,"
                "mass_local_linf,mass_normalized_imbalance,"
                "nonfinite_mass_faces,nonfinite_velocities,"
                "high_order_fraction,face_momentum_x_relative,"
                "face_momentum_y_relative,face_momentum_z_relative,"
                "face_momentum_x_linf,face_momentum_y_linf,"
                "face_momentum_z_linf,convection_scheme,"
                "momentum_x_rhs_relative,momentum_y_rhs_relative,"
                "momentum_z_rhs_relative,momentum_component_max,"
                "pressure_residual,velocity_change_inf,pressure_change_inf\n"
                "1,1e-11,1e-12,2e-11,1e-13,0,0,1,1e-11,1e-11,0,1e-11,1e-11,0,"
                f"{scheme},1e-11,1e-11,0,1e-11,1e-12,1e-12,1e-12\n"
            )
    result = analyze(tmp_path / "families", tmp_path / "runs", Path(sys.executable))
    assert result["status"] == expected
    if expected == "BLOCKED":
        assert any(
            "incompatible physical setup: convection_scheme" in reason
            for reason in result["reasons"]
        )
    assert result["causal_anisotropy_isolated"] is False
    assert result["comparison_scope"] == (
        "TWO_MESH_DISTRIBUTION_FAMILIES_AND_POLYGONAL_WALL_APPROXIMATIONS"
    )
    for name, scheme in zip(names, schemes):
        assert result["families"][name]["status"] == "GCI_CANDIDATE"
        assert result["families"][name]["physical_setup"]["convection_scheme"] == scheme
        for component in ("cd", "cdp", "cdv"):
            assert (
                result["families"][name]["gci"][component]["status"]
                == "asymptotic_candidate"
            )


def test_face_flux_restart_pair_is_fail_closed(tmp_path):
    mesh = tmp_path / "mesh.h5"
    mesh.write_bytes(b"placeholder")
    flux = tmp_path / "massflux.csv"
    flux.write_text("placeholder")
    with pytest.raises(ValueError, match="paired same-mesh cell state"):
        run_case(
            mesh,
            tmp_path / "result",
            Path(sys.executable),
            10,
            10,
            initial_face_flux=flux,
        )
    state = tmp_path / "state.csv"
    state.write_text("placeholder")
    with pytest.raises(FileNotFoundError, match="face-flux restart file"):
        run_case(
            mesh,
            tmp_path / "result",
            Path(sys.executable),
            10,
            10,
            initial_state=state,
            initial_face_flux=tmp_path / "nonexistent.csv",
        )


def test_pressure_bc_policy_must_be_explicit_and_flux_certified(tmp_path):
    """Direct run-metadata preflight protects against silent mixing at GCI entry."""
    name = "quasi_isotropic"
    folder = tmp_path / "families" / name
    folder.mkdir(parents=True)
    rows = quality_rows()
    (folder / "quality.json").write_text(json.dumps(rows))
    for n in (64, 128, 256):
        location = tmp_path / "runs" / name / f"n{n}"
        location.mkdir(parents=True)
        (location / "run_summary.json").write_text(
            json.dumps(
                {
                    "mesh": str((folder / f"naca0012_{name}_n{n}.h5").resolve()),
                    "scheme": "TVD_VANLEER",
                    "algorithm": "SIMPLE",
                    "linear_tolerance": 1e-12,
                    "wall_traction": "FACE_WLS",
                    "nonlinear_tolerance": 1e-8,
                    "exit_code": 0,
                    "kinematic_face_flux_certified": True,
                    "pressure_boundary_policy": "ZERO_NORMAL_GRADIENT",
                }
            )
        )
    # No force/history fixtures needed: a different boundary contract should
    # be rejected before any physical certificate is inspected.
    target = tmp_path / "runs" / name / "n128" / "run_summary.json"
    base = json.loads(target.read_text())
    for mutation in (
        {"pressure_boundary_policy": "DIRICHLET_ZERO"},
        {"pressure_boundary_policy": None},
        {"kinematic_face_flux_certified": False},
        {"algorithm": "COUPLED"},
        {"algorithm": None},
    ):
        altered = {**base, **mutation}
        target.write_text(json.dumps(altered))
        result = GCI.inspect_family(
            name, folder, tmp_path / "runs", Path(sys.executable)
        )
        assert result["status"] == "BLOCKED"
        assert any("N128:" in reason for reason in result["reasons"]) or (
            "mixed farfield pressure policies across N64/N128/N256" in result["reasons"]
        )


def test_two_grid_provisional_cannot_mix_vanleer_and_muscl(tmp_path):
    """Even assumed-order, explicitly non-certified output must not mix schemes."""
    name = "quasi_isotropic"
    folder = tmp_path / "families" / name
    folder.mkdir(parents=True)
    (folder / "quality.json").write_text(json.dumps(quality_rows()))
    for n, scheme in ((64, "TVD_VANLEER"), (128, "MUSCL_CELL_LIMITED")):
        run = tmp_path / "runs" / name / f"n{n}"
        run.mkdir(parents=True)
        (run / "run_summary.json").write_text(json.dumps({
            "mesh": str((folder / f"naca0012_{name}_n{n}.h5").resolve()),
            "scheme": scheme, "algorithm": "SIMPLE",
            "linear_tolerance": 1e-12, "wall_traction": "FACE_WLS",
            "nonlinear_tolerance": 1e-8, "exit_code": 0,
            "pressure_boundary_policy": "ZERO_NORMAL_GRADIENT",
            "kinematic_face_flux_certified": True,
        }))
        (run / "force.csv").write_text(
            "case,h,converged,res_mom,res_mass,cd,cdp,cdv,cl,cm\n"
            f"N{n},0.01,yes,1e-11,1e-11,0.12,0.04,0.08,0,0\n")
        write_synthetic_raw_wall_faces(run / "raw_wall_faces.csv", 0.04, 0.08)
        (run / "history.csv").write_text(
            "iteration,momentum_rhs_relative,continuity_normalized,"
            "mass_local_linf,mass_normalized_imbalance,nonfinite_mass_faces,"
            "nonfinite_velocities,high_order_fraction,face_momentum_x_relative,"
            "face_momentum_y_relative,face_momentum_z_relative,face_momentum_x_linf,"
            "face_momentum_y_linf,face_momentum_z_linf,convection_scheme,"
            "momentum_x_rhs_relative,momentum_y_rhs_relative,momentum_z_rhs_relative,"
            "momentum_component_max,pressure_residual,velocity_change_inf,"
            "pressure_change_inf\n"
            "1,1e-11,1e-12,2e-11,1e-13,0,0,1,1e-11,1e-11,0,1e-11,1e-11,0,"
            f"{scheme},1e-11,1e-11,0,1e-11,1e-12,1e-12,1e-12\n")
    result = GCI.inspect_family(name, folder, tmp_path / "runs", Path(sys.executable))
    assert result["status"] == "BLOCKED"
    assert "two_grid_assumed_order_only" not in result
    assert any("mixed final convection schemes" in r for r in result["reasons"])


def test_two_grid_provisional_cannot_mix_vanleer_and_muscl(tmp_path):
    """Even assumed-order, explicitly non-certified output must not mix schemes."""
    name = "quasi_isotropic"
    folder = tmp_path / "families" / name
    folder.mkdir(parents=True)
    (folder / "quality.json").write_text(json.dumps(quality_rows()))
    for n, scheme in ((64, "TVD_VANLEER"), (128, "MUSCL_CELL_LIMITED")):
        run = tmp_path / "runs" / name / f"n{n}"
        run.mkdir(parents=True)
        (run / "run_summary.json").write_text(json.dumps({
            "mesh": str((folder / f"naca0012_{name}_n{n}.h5").resolve()),
            "scheme": scheme, "algorithm": "SIMPLE",
            "linear_tolerance": 1e-12, "wall_traction": "FACE_WLS",
            "nonlinear_tolerance": 1e-8, "exit_code": 0,
            "pressure_boundary_policy": "ZERO_NORMAL_GRADIENT",
            "kinematic_face_flux_certified": True,
        }))
        (run / "force.csv").write_text(
            "case,h,converged,res_mom,res_mass,cd,cdp,cdv,cl,cm\n"
            f"N{n},0.01,yes,1e-11,1e-11,0.12,0.04,0.08,0,0\n")
        write_synthetic_raw_wall_faces(run / "raw_wall_faces.csv", 0.04, 0.08)
        (run / "history.csv").write_text(
            "iteration,momentum_rhs_relative,continuity_normalized,"
            "mass_local_linf,mass_normalized_imbalance,nonfinite_mass_faces,"
            "nonfinite_velocities,high_order_fraction,face_momentum_x_relative,"
            "face_momentum_y_relative,face_momentum_z_relative,face_momentum_x_linf,"
            "face_momentum_y_linf,face_momentum_z_linf,convection_scheme,"
            "momentum_x_rhs_relative,momentum_y_rhs_relative,momentum_z_rhs_relative,"
            "momentum_component_max,pressure_residual,velocity_change_inf,"
            "pressure_change_inf\n"
            "1,1e-11,1e-12,2e-11,1e-13,0,0,1,1e-11,1e-11,0,1e-11,1e-11,0,"
            f"{scheme},1e-11,1e-11,0,1e-11,1e-12,1e-12,1e-12\n")
    result = GCI.inspect_family(name, folder, tmp_path / "runs", Path(sys.executable))
    assert result["status"] == "BLOCKED"
    assert "two_grid_assumed_order_only" not in result
    assert any("mixed final convection schemes" in r for r in result["reasons"])


def test_two_grid_provisional_cannot_mix_vanleer_and_muscl(tmp_path):
    """Even assumed-order, explicitly non-certified output must not mix schemes."""
    name = "quasi_isotropic"
    folder = tmp_path / "families" / name
    folder.mkdir(parents=True)
    (folder / "quality.json").write_text(json.dumps(quality_rows()))
    for n, scheme in ((64, "TVD_VANLEER"), (128, "MUSCL_CELL_LIMITED")):
        run = tmp_path / "runs" / name / f"n{n}"
        run.mkdir(parents=True)
        (run / "run_summary.json").write_text(
            json.dumps(
                {
                    "mesh": str((folder / f"naca0012_{name}_n{n}.h5").resolve()),
                    "scheme": scheme,
                    "algorithm": "SIMPLE",
                    "linear_tolerance": 1e-12,
                    "wall_traction": "FACE_WLS",
                    "nonlinear_tolerance": 1e-8,
                    "exit_code": 0,
                    "pressure_boundary_policy": "ZERO_NORMAL_GRADIENT",
                    "kinematic_face_flux_certified": True,
                }
            )
        )
        (run / "force.csv").write_text(
            "case,h,converged,res_mom,res_mass,cd,cdp,cdv,cl,cm\n"
            f"N{n},0.01,yes,1e-11,1e-11,0.12,0.04,0.08,0,0\n"
        )
        write_synthetic_raw_wall_faces(run / "raw_wall_faces.csv", 0.04, 0.08)
        (run / "history.csv").write_text(
            "iteration,momentum_rhs_relative,continuity_normalized,"
            "mass_local_linf,mass_normalized_imbalance,nonfinite_mass_faces,"
            "nonfinite_velocities,high_order_fraction,face_momentum_x_relative,"
            "face_momentum_y_relative,face_momentum_z_relative,face_momentum_x_linf,"
            "face_momentum_y_linf,face_momentum_z_linf,convection_scheme,"
            "momentum_x_rhs_relative,momentum_y_rhs_relative,momentum_z_rhs_relative,"
            "momentum_component_max,pressure_residual,velocity_change_inf,"
            "pressure_change_inf\n"
            "1,1e-11,1e-12,2e-11,1e-13,0,0,1,1e-11,1e-11,0,1e-11,1e-11,0,"
            f"{scheme},1e-11,1e-11,0,1e-11,1e-12,1e-12,1e-12\n"
        )
    result = GCI.inspect_family(name, folder, tmp_path / "runs", Path(sys.executable))
    assert result["status"] == "BLOCKED"
    assert "two_grid_assumed_order_only" not in result
    assert any("mixed final convection schemes" in r for r in result["reasons"])
