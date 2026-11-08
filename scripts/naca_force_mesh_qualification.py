"""Fail-closed three-grid NACA force qualification with actual mesh evidence.

The solver must first produce a converged force certificate for each mesh.
The independent nonlinear history must agree with the final certificate,
including a separate mass-conservation gate. The mesh-quality report must
prove comparable tangential AND near-wall refinement before an ordinary GCI.

No preliminary, failed, missing or stale solution is ever promoted to a GCI.
This is grid VERIFICATION, not validation against experimental forces.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

from naca_mesh_similarity_gate import verify_mesh_family
from naca_surface_force_oracle import integrate as integrate_raw_faces
from naca_surface_force_oracle import read_rows as read_raw_faces
from qualify_naca_grid import REQUIRED, qualify

LEVELS = (64, 128, 256)
RELATIVE_GATE = 1.0e-8


def read_one_certificate(path: Path, level: int) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if set(reader.fieldnames or ()) != set(REQUIRED):
            raise ValueError(f"N{level}: force certificate has wrong CSV fields")
        rows = list(reader)
    if len(rows) != 1:
        raise ValueError(f"N{level}: exactly one force certificate row required")
    row = rows[0]
    if row["case"] != f"N{level}":
        raise ValueError(f"N{level}: force certificate uses wrong mesh label")
    if row["converged"].lower() not in ("yes", "true", "1"):
        raise ValueError(f"N{level}: nonlinear solve not converged")
    for key in ("h", "res_mom", "res_mass", "cd", "cdp", "cdv", "cl", "cm"):
        if not math.isfinite(float(row[key])):
            raise ValueError(f"N{level}: nonfinite {key}")
    if float(row["res_mom"]) < 0 or float(row["res_mass"]) < 0:
        raise ValueError(f"N{level}: negative residual norm in force certificate")
    if abs(float(row["h"]) - 1.0 / level) > 1e-12:
        raise ValueError(f"N{level}: wrong characteristic tangential spacing")
    if max(float(row["res_mom"]), float(row["res_mass"])) > RELATIVE_GATE:
        raise ValueError(f"N{level}: nonlinear momentum or mass gate not met")
    if abs(float(row["cd"]) - float(row["cdp"]) - float(row["cdv"])) > 1e-8:
        raise ValueError(f"N{level}: pressure and viscous forces do not sum to Cd")
    if abs(float(row["cl"])) > 2e-2 or float(row["cdp"]) <= 0.0:
        raise ValueError(f"N{level}: sign or symmetry physical gate failed")
    return row


def verify_physical_history(
    path: Path, level: int, certificate: dict[str, str]
) -> dict:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        # Matrix reassembly and conservative face flux re-evaluation are two
        # different residuals. Qualification needs both on the accepted state.
        independent = {
            f"face_momentum_{axis}_{field}"
            for axis in ("x", "y", "z")
            for field in ("relative", "linf")
        }
        needed = {
            "iteration",
            "momentum_rhs_relative",
            "continuity_normalized",
            "mass_local_linf",
            "mass_normalized_imbalance",
            "nonfinite_mass_faces",
            "nonfinite_velocities",
            "high_order_fraction",
            "convection_scheme",
            # Independent checks that the solver itself requires. Refuse an
            # external/corrupted certificate masking a small-component error.
            "momentum_x_rhs_relative",
            "momentum_y_rhs_relative",
            "momentum_z_rhs_relative",
            "momentum_component_max",
            "pressure_residual",
            "velocity_change_inf",
            "pressure_change_inf",
        } | independent
        if not needed.issubset(reader.fieldnames or ()):
            raise ValueError(f"N{level}: incomplete physical history schema")
        count = 0
        last = None
        run_scheme = None
        for row in reader:
            count += 1
            iteration = int(row["iteration"])
            if iteration != count:
                raise ValueError(
                    f"N{level}: incomplete or nonconsecutive history iteration"
                )
            for key in needed - {"iteration", "convection_scheme"}:
                if not math.isfinite(float(row[key])):
                    raise ValueError(f"N{level}: nonfinite history {key}")
            for key in {
                "momentum_rhs_relative",
                "continuity_normalized",
                "mass_local_linf",
                "mass_normalized_imbalance",
                "momentum_x_rhs_relative",
                "momentum_y_rhs_relative",
                "momentum_z_rhs_relative",
                "momentum_component_max",
                "pressure_residual",
                "velocity_change_inf",
                "pressure_change_inf",
            } | independent:
                if float(row[key]) < 0:
                    raise ValueError(f"N{level}: negative physical residual norm {key}")
            if not 0.0 <= float(row["high_order_fraction"]) <= 1.0:
                raise ValueError(f"N{level}: invalid convection fraction")
            components = [float(row[f"momentum_{axis}_rhs_relative"]) for axis in "xyz"]
            if not math.isclose(
                float(row["momentum_component_max"]),
                max(components),
                abs_tol=1e-11,
                rel_tol=1e-6,
            ):
                raise ValueError(f"N{level}: inconsistent momentum component maximum")
            # A startup continuation may legitimately use UPWIND. Only the
            # accepted FINAL state must be fully second-order bounded TVD;
            # never reject a valid startup solely for its early scheme.
            if row["convection_scheme"] not in (
                "UPWIND",
                "TVD",
                "TVD_VANLEER",
                "TVD_VANALBADA",
                "MUSCL_CELL_LIMITED",
            ):
                raise ValueError(f"N{level}: unknown convection scheme")
            if run_scheme is None:
                run_scheme = row["convection_scheme"]
            elif row["convection_scheme"] != run_scheme:
                raise ValueError(
                    f"N{level}: nonlinear history mixes convection schemes"
                )
            if (
                int(row["nonfinite_mass_faces"]) != 0
                or int(row["nonfinite_velocities"]) != 0
            ):
                raise ValueError(f"N{level}: nonfinite conservation or velocity")
            last = row
    if count == 0 or last is None:
        raise ValueError(f"N{level}: missing physical history")
    momentum = float(last["momentum_rhs_relative"])
    mass_local = float(last["mass_local_linf"])
    mass_imbalance = float(last["mass_normalized_imbalance"])
    continuity = float(last["continuity_normalized"])
    # Must match the C++ force certificate writer exactly.
    mass = max(continuity, mass_imbalance, mass_local)
    if momentum > RELATIVE_GATE or mass > RELATIVE_GATE:
        raise ValueError(f"N{level}: nonlinear history is not physically converged")
    if any(
        float(last[key]) > RELATIVE_GATE
        for key in (
            "momentum_component_max",
            "pressure_residual",
            "velocity_change_inf",
            "pressure_change_inf",
        )
    ):
        raise ValueError(
            f"N{level}: transverse momentum/linear/change convergence gate not met"
        )
    if float(last["high_order_fraction"]) != 1.0:
        raise ValueError(f"N{level}: final convection operator is not full order")
    if last["convection_scheme"] not in (
        "TVD",
        "TVD_VANLEER",
        "TVD_VANALBADA",
        "MUSCL_CELL_LIMITED",
    ):
        raise ValueError(
            f"N{level}: final qualification requires bounded "
            f"second-order TVD or MUSCL convection"
        )
    face_residuals = {
        axis: float(last[f"face_momentum_{axis}_relative"]) for axis in ("x", "y", "z")
    }
    if any(value > RELATIVE_GATE for value in face_residuals.values()):
        raise ValueError(
            f"N{level}: independent local face momentum residual is not converged"
        )
    for key, value in (("res_mom", momentum), ("res_mass", mass)):
        claimed = float(certificate[key])
        if abs(claimed - value) > 1e-13 + 1e-5 * max(abs(claimed), abs(value)):
            raise ValueError(f"N{level}: stale or mismatched force certificate ({key})")
    return {
        "iterations": count,
        "final_momentum_rhs_relative": momentum,
        "final_momentum_component_max": float(last["momentum_component_max"]),
        "final_continuity_normalized": continuity,
        "final_mass_local_linf": mass_local,
        "final_mass_normalized_imbalance": mass_imbalance,
        "final_certificate_mass_residual": mass,
        "final_independent_face_momentum_relative": face_residuals,
        "convection_scheme": run_scheme,
        "experimental_reconstruction": run_scheme == "MUSCL_CELL_LIMITED",
    }


def certify(
    mesh_quality: list[dict],
    certs: dict[int, Path],
    histories: dict[int, Path],
    raw_faces: dict[int, Path] | None = None,
) -> dict:
    mesh_report = verify_mesh_family(mesh_quality)
    by_n = {int(row["n_surface"]): row for row in mesh_quality}
    if set(by_n) != set(LEVELS):
        raise ValueError("mesh quality data must describe exactly N64, N128, N256")
    report: dict = {
        "status": "BLOCKED",
        "mesh_family": mesh_report,
        "cases": [],
        "reasons": [],
        "gci": {"status": "NOT_QUALIFIED", "validation_reference": "NOT_PERFORMED"},
    }
    rows: list[dict[str, str]] = []
    schemes: set[str] = set()
    if raw_faces is None:
        report["reasons"].append(
            "missing independent raw wall-face traction files for N64/N128/N256"
        )
    for level in LEVELS:
        try:
            certificate = read_one_certificate(certs[level], level)
            physical = verify_physical_history(histories[level], level, certificate)
            schemes.add(str(physical["convection_scheme"]))
            # An independent script must recompute the face-area/normal,
            # pressure/viscous and moment integrals from *raw* traction.
            # This still does not validate the underlying wall stress model.
            if raw_faces is None:
                raise ValueError(f"N{level}: missing independent raw force integration")
            face_result = integrate_raw_faces(
                read_raw_faces(raw_faces[level]), certificate
            )
            rows.append(certificate)
            report["cases"].append(
                {
                    "N": level,
                    "status": "CONVERGED",
                    "Cd": float(certificate["cd"]),
                    "Cd_pressure": float(certificate["cdp"]),
                    "Cd_viscous": float(certificate["cdv"]),
                    "independent_face_integration": face_result,
                    **physical,
                }
            )
        except (OSError, ValueError, KeyError, TypeError, OverflowError) as exc:
            report["cases"].append({"N": level, "status": "NOT_QUALIFIED"})
            report["reasons"].append(str(exc))
    # Richardson/GCI estimates are meaningless when each level discretizes a
    # different convection operator, even if all three solved to tight norms.
    if len(schemes) > 1:
        report["reasons"].append(
            "three-grid Richardson study mixes different convection schemes"
        )
    if len(schemes) == 1:
        report["convection_scheme"] = next(iter(schemes))
        report["experimental_reconstruction"] = (
            report["convection_scheme"] == "MUSCL_CELL_LIMITED"
        )
    if mesh_report["status"] != "COMPATIBLE_WITH_SINGLE_RATIO_GCI":
        report["reasons"].append(
            "near-wall/tangential grid refinement is not geometrically similar"
        )
    if report["reasons"]:
        return report
    try:
        gci = qualify(rows, RELATIVE_GATE)
        # Total drag can look monotone even when opposing pressure and skin-
        # friction errors cancel. Each constituent must independently exhibit
        # a usable three-grid Richardson trend before force GCI promotion.
        for component in ("cd", "cdp", "cdv"):
            if gci["quantities"][component]["status"] != "asymptotic_candidate":
                raise ValueError(
                    f"{component}: drag component lacks an asymptotic Richardson signature"
                )
        report["gci"] = {
            "status": "RICHARDSON_CANDIDATE_ONLY",
            "result": gci,
            "validation_reference": "NOT_PERFORMED",
        }
        report["status"] = "GCI_CANDIDATE"
    except (ValueError, KeyError, ZeroDivisionError) as exc:
        report["reasons"].append(str(exc))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-quality", required=True, type=Path)
    for level in LEVELS:
        parser.add_argument(f"--cert-{level}", required=True, type=Path)
        parser.add_argument(f"--history-{level}", required=True, type=Path)
        parser.add_argument(f"--raw-faces-{level}", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        quality = json.loads(args.mesh_quality.read_text(encoding="utf-8"))
        certs = {level: getattr(args, f"cert_{level}") for level in LEVELS}
        histories = {level: getattr(args, f"history_{level}") for level in LEVELS}
        raw_faces = {level: getattr(args, f"raw_faces_{level}") for level in LEVELS}
        result = certify(quality, certs, histories, raw_faces)
        # Replacement only after a fully written JSON, never a half-written
        # or stale SUCCESS from an older campaign.
        temporary = args.output.with_name(args.output.name + ".tmp")
        temporary.write_text(
            json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
        )
        os.replace(temporary, args.output)
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        args.output.unlink(missing_ok=True)
        print(f"NACA FORCE MESH QUALIFICATION FAILED: {exc}")
        return 2
    print(json.dumps(result, sort_keys=True, indent=2))
    if result["status"] != "GCI_CANDIDATE":
        print("NACA GRID GCI: BLOCKED (no physically verified 3-grid result)")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
