"""Independent NACA momentum bottleneck investigation (diagnostic; never a GCI gate).

Inputs MUST describe the same accepted iteration. The script recomputes cell
mass divergence solely from the saved authoritative face flux, checks against
the independently recorded history, and inspects the bounded-versus-physical
momentum balance. A centreline outflow face is used to quantify (not to prove
causation from) the imposed wake-to-freestream jump.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

TOL = 1.0e-8


def read_csv(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = csv.DictReader(stream)
        if rows.fieldnames is None or any(key not in rows.fieldnames for key in required):
            raise ValueError(f"{path}: missing required CSV columns")
        result = list(rows)
    if not result:
        raise ValueError(f"{path}: empty CSV")
    return result


def number(row: dict[str, str], key: str) -> float:
    v = float(row[key])
    if not math.isfinite(v):
        raise ValueError(f"nonfinite {key}: {v}")
    return v


def iteration_rows(history: list[dict[str, str]], budget: list[dict[str, str]]):
    iterations = [int(number(r, "iteration")) for r in history]
    if sorted(iterations) != iterations or len(set(iterations)) != len(iterations):
        raise ValueError("history iterations are not strictly increasing")
    b_iters = {int(number(r, "iteration")) for r in budget}
    if len(b_iters) != 1 or next(iter(b_iters)) != iterations[-1]:
        raise ValueError("refuse to mix different accepted iterations: history/budget")
    return history[-1]


def mass_flux_audit(cells: list[dict[str, str]],
                    faces: list[dict[str, str]],
                    density: float = 1.0,
                    u_inf: float = 1.0,
                    outlet_x_min: float = 19.0) -> dict[str, float]:
    if not (density > 0 and u_inf > 0):
        raise ValueError("density and freestream velocity must be positive")
    by_cell = {int(number(c, "cell")): c for c in cells}
    if len(by_cell) != len(cells):
        raise ValueError("duplicated cell ids")
    mass = {idx: 0.0 for idx in by_cell}
    total_abs_flux = 0.0
    best = None
    best_score = math.inf
    for f in faces:
        owner = int(number(f, "owner"))
        neighbour = int(number(f, "neighbour"))
        flux = number(f, "mass_flux")
        if owner not in mass or (neighbour >= 0 and neighbour not in mass):
            raise ValueError("invalid face connectivity")
        mass[owner] += flux
        if neighbour >= 0:
            mass[neighbour] -= flux
        total_abs_flux += abs(flux)
        sx = number(f, "sx")
        sy = number(f, "sy")
        sz = number(f, "sz")
        area = math.sqrt(sx*sx + sy*sy + sz*sz)
        if (neighbour < 0 and number(f, "cx") > outlet_x_min and
            area > 0 and sx / area > 0.90):
            distance_to_centre = abs(number(f, "cy"))
            if distance_to_centre < best_score:
                best = (f, owner)
                best_score = distance_to_centre
    if best is None:
        raise ValueError("no downstream boundary face in requested region")
    if total_abs_flux <= 0:
        raise ValueError("zero flux norm")
    face, owner = best
    state = by_cell[owner]
    actual_flux = number(face, "mass_flux")
    prescribed_flux = density*u_inf*number(face, "sx")
    extrapolated_flux = density*(
        number(state, "ux")*number(face, "sx") +
        number(state, "uy")*number(face, "sy"))
    return {
        "cell_count": len(by_cell), "face_count": len(faces),
        "continuity_l1": sum(abs(x) for x in mass.values()),
        "continuity_linf": max(abs(x) for x in mass.values()),
        "continuity_l1_relative": sum(abs(x) for x in mass.values()) / total_abs_flux,
        "sum_abs_mass_flux": total_abs_flux,
        "continuity_worst_cell": max(mass, key=lambda c: abs(mass[c])),
        "centreline_outlet_face": int(number(face, "face")),
        "centreline_outlet_owner": owner,
        "centreline_outlet_x": number(face, "cx"),
        "centreline_outlet_ux": number(state, "ux"),
        "centreline_outlet_uy": number(state, "uy"),
        "prescribed_face_flux": prescribed_flux,
        "authoritative_face_flux": actual_flux,
        "cell_velocity_extrapolated_flux": extrapolated_flux,
        "boundary_flux_contract_error": actual_flux - prescribed_flux,
        "wake_to_prescribed_flux_gap": actual_flux - extrapolated_flux,
        "wake_to_prescribed_relative_gap":
            (actual_flux - extrapolated_flux) / max(abs(prescribed_flux), 1e-300),
    }


def worst_budget(budget: list[dict[str, str]], component: int,
                 cells: list[dict[str, str]]) -> dict[str, float]:
    by_cell = {int(number(c, "cell")): c for c in cells}
    matching = [b for b in budget if int(number(b, "component")) == component]
    if not matching:
        raise ValueError(f"missing budget component {component}")
    b = max(matching, key=lambda r: abs(number(r, "physical_residual")))
    phys = number(b, "physical_residual")
    bounded = number(b, "bounded_operator_residual")
    correction = number(b, "bounded_divergence")
    mass_div = number(b, "mass_divergence")
    if not math.isclose(bounded, phys + correction, abs_tol=1e-12, rel_tol=1e-9):
        raise ValueError("physical/bounded budget algebra is inconsistent")
    cell_id = int(number(b, "cell"))
    velocity_key = ("ux", "uy", "uz")[component]
    if velocity_key in by_cell[cell_id]:
        velocity = number(by_cell[cell_id], velocity_key)
        if not math.isclose(correction, -velocity*mass_div,
                            abs_tol=1e-12, rel_tol=1e-6):
            raise ValueError("bounded divergence differs from -U*div(mass_flux)")
    conv = number(b, "convection")
    pressure = number(b, "pressure")
    visc = number(b, "viscous") + number(b, "transpose_viscous")
    return {
        "cell": cell_id,
        "x": number(b, "x"), "y": number(b, "y"),
        "convection": conv, "pressure": pressure, "viscosity": visc,
        "mass_divergence": mass_div,
        "bounded_correction": correction,
        "physical_residual": phys,
        "bounded_operator_residual": bounded,
        "fraction_physical_from_mass":
            abs(correction) / max(abs(phys), 1e-300),
        "relative_uncancelled_balance":
            abs(phys) / max(abs(conv) + abs(visc) + abs(pressure), 1e-300),
    }


def log_decay(history: list[dict[str, str]], key: str,
              last_n: int = 160) -> dict[str, float]:
    rows = history[-min(last_n, len(history)):]
    xx = [number(r, "iteration") for r in rows]
    yy = [math.log(number(r, key)) for r in rows]
    if len(xx) < 8:
        raise ValueError("insufficient history for decay estimate")
    x0 = sum(xx)/len(xx)
    y0 = sum(yy)/len(yy)
    var = sum((x-x0)**2 for x in xx)
    slope = sum((x-x0)*(y-y0) for x,y in zip(xx,yy)) / var
    predicted = [y0 + slope*(x-x0) for x in xx]
    ss_tot = sum((y-y0)**2 for y in yy)
    r2 = 1 - sum((y-p)**2 for y,p in zip(yy,predicted)) / ss_tot if ss_tot else 0.0
    final_value = number(rows[-1], key)
    remaining = (
        math.log(TOL/final_value)/slope
        if slope < 0 and final_value > TOL and r2 >= .90 else None
    )
    return {"first": number(rows[0], key), "last": final_value,
            "log_slope_per_iteration": slope, "fit_r2": r2,
            "hypothetical_additional_iterations":
                math.ceil(remaining) if remaining is not None else None}



def checkpoint_metadata(path: Path, state_path: Path, flux_path: Path,
                        cell_count: int, face_count: int) -> dict[str, int]:
    """Cross-check an accepted V3 manifest against the actual snapshot files."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 2 or lines[0].strip() != "CFDX_NACA_PAIRED_V3":
        raise ValueError("checkpoint is not a paired V3 snapshot")
    parts = lines[1].split(",")
    if len(parts) != 10:
        raise ValueError("invalid V3 manifest fields")
    absolute_iteration, cells, faces = [int(x) for x in parts[:3]]
    if (absolute_iteration < 1 or cells != cell_count or faces != face_count or
        parts[5] != state_path.name or parts[6] != flux_path.name):
        raise ValueError("checkpoint manifest does not match state/flux pair")
    return {"accepted_iteration_absolute": absolute_iteration,
            "cell_count": cells, "face_count": faces}


def spatial_residual_signature(budget: list[dict[str, str]], component: int,
                               top_n: int = 10) -> dict[str, float]:
    """Locate the largest physical residuals without claiming causation."""
    candidates = [b for b in budget if int(number(b,"component")) == component]
    if not candidates:
        raise ValueError("no matching momentum components")
    candidates.sort(key=lambda row: abs(number(row,"physical_residual")),reverse=True)
    top = candidates[:top_n]
    mirror_pair = False
    if component == 1 and len(top) >= 2:
        a,b = top[:2]
        mirror_pair = (abs(number(a,"x")-number(b,"x")) < 1e-7 and
                       abs(number(a,"y")+number(b,"y")) < 1e-7 and
                       abs(number(a,"physical_residual")+
                           number(b,"physical_residual")) < 1e-12)
    return {
        "top_count": len(top),
        "count_at_x_ge_15": sum(number(row,"x") >= 15 for row in top),
        "count_at_x_ge_19_5": sum(number(row,"x") >= 19.5 for row in top),
        "worst_two_mirror_symmetric": mirror_pair,
        "worst_x": number(top[0],"x"),
        "worst_y": number(top[0],"y"),
    }

def diagnose(history_path: Path, budget_path: Path, state_path: Path,
             flux_path: Path, long_history_path: Path | None = None,
             manifest_path: Path | None = None) -> dict:
    history = read_csv(history_path, ("iteration", "continuity_normalized",
                                     "momentum_y_rhs_relative", "momentum_reference_relative"))
    budget = read_csv(budget_path, ("iteration","component","cell",
                                   "bounded_divergence","physical_residual",
                                   "bounded_operator_residual"))
    final = iteration_rows(history, budget)
    cells = read_csv(state_path, ("cell","ux","uy"))
    faces = read_csv(flux_path, ("face","owner","neighbour",
                                 "cx","cy","sx","sy","sz","mass_flux"))
    audit = mass_flux_audit(cells, faces)
    if not math.isclose(audit["continuity_l1_relative"], number(final,"continuity_normalized"),
                        rel_tol=2e-5, abs_tol=1e-12):
        raise ValueError("saved face flux does not match history continuity")
    if "mass_local_linf" in final and not math.isclose(
            audit["continuity_linf"], number(final,"mass_local_linf"),
            rel_tol=2e-5, abs_tol=1e-12):
        raise ValueError("saved face flux does not match worst-cell mass imbalance")
    result = {"same_accepted_iteration": int(number(final,"iteration")),
              "history": {k: number(final,k) for k in (
                  "momentum_rhs_relative","momentum_y_rhs_relative",
                  "momentum_reference_relative","face_momentum_y_relative",
                  "continuity_normalized","velocity_change_inf",
                  "pressure_change_inf","pressure_residual")},
              "mass_and_outlet": audit,
              "worst_x": worst_budget(budget, 0, cells),
              "worst_y": worst_budget(budget, 1, cells),
              "spatial_x": spatial_residual_signature(budget, 0),
              "spatial_y": spatial_residual_signature(budget, 1)}
    if manifest_path:
        result["manifest"] = checkpoint_metadata(manifest_path, state_path,
                                                  flux_path, len(cells), len(faces))
    if long_history_path:
        long = read_csv(long_history_path, ("iteration","momentum_y_rhs_relative"))
        result["long_run_rates"] = {
            k: log_decay(long, k) for k in (
                "momentum_y_rhs_relative","momentum_reference_relative",
                "face_momentum_y_relative","velocity_change_inf",
                "pressure_change_inf")}
        result["long_run_iteration_is_not_same_as_snapshot"] = True
    result["physically_qualified"] = False  # This diagnostic NEVER issues a certificate.
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ("history", "budget", "state", "flux"):
        ap.add_argument("--"+key, type=Path, required=True)
    ap.add_argument("--long-history", type=Path)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    result = diagnose(args.history, args.budget, args.state, args.flux,
                      args.long_history, args.manifest)
    output = json.dumps(result, indent=2, allow_nan=False, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()

