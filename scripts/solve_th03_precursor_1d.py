#!/usr/bin/env python3
"""Experimental 1-D *fully developed channel* SST-2003 precursor for TH-03.

Solves finite-volume balances for U, k and omega at fixed bulk flow.
This is deliberately NOT a measured Driver--Seegmiller inflow profile, nor
an x-marching growing boundary layer. Its discrete equations, corrections,
iteration changes, and bulk flow must all pass before exporting the state.

Smooth-wall FV boundary conditions use U=k=0 and the original CFDX omega
face Dirichlet convention 6 nu/(beta1*d_1**2). It never silently substitutes
the unrelated NASA first-centre factor 60 prescription.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded

from scripts.generate_th03_openfoam_blockmesh import NU, VELOCITY, H


@dataclass(frozen=True)
class PrecursorControls:
    n_cells: int = 192
    stretch: float = 4.0
    max_iterations: int = 2500
    tolerance: float = 2e-7
    damping: float = 0.25
    bulk_velocity: float = VELOCITY
    height: float = 8.0 * H
    beta_star: float = 0.09
    beta1: float = 0.075
    beta2: float = 0.0828
    gamma1: float = 5.0 / 9.0
    gamma2: float = 0.44
    a1: float = 0.31
    sigma_k1: float = 0.85
    sigma_k2: float = 1.0
    sigma_w1: float = 0.5
    sigma_w2: float = 0.856
    production_limiter: float = 10.0
    k_floor: float = 1.0e-10
    omega_floor: float = 1.0e-10


def grid(c: PrecursorControls):
    if (
        c.n_cells < 8
        or c.n_cells > 4096
        or not (0.0 <= c.stretch < 10.0)
        or not (0.0 < c.damping <= 1.0)
        or not (0.0 < c.tolerance <= 1e-3)
        or c.max_iterations < 1
        or not (c.bulk_velocity > 0.0 and c.height > 0.0)
    ):
        raise ValueError("invalid 1D precursor controls")
    xi = np.linspace(0.0, 1.0, c.n_cells + 1)
    if c.stretch > 1.0e-12:
        edges = c.height * (
            0.5 + 0.5 * np.tanh(c.stretch * (2.0 * xi - 1.0)) / np.tanh(c.stretch)
        )
    else:
        edges = c.height * xi
    y = 0.5 * (edges[:-1] + edges[1:])
    return edges, y


def fv_diffusion(edges, y, gamma, sink, source, bc0, bc1, *, wall_gamma=None):
    """Return SPD M-matrix solution for -div(G grad f)+S_p*f=S_u.

    Coefficients are integrated over each finite volume; the residual uses
    the SAME physical FV operator including the exact boundary source.
    """
    n = len(y)
    gamma = np.asarray(gamma, float)
    sink = np.asarray(sink, float)
    source = np.asarray(source, float)
    if not (
        gamma.shape == sink.shape == source.shape == (n,)
        and np.all(np.isfinite(gamma))
        and np.all(gamma > 0)
        and np.all(np.isfinite(sink))
        and np.all(sink >= 0)
        and np.all(np.isfinite(source))
        and np.isfinite(bc0)
        and np.isfinite(bc1)
    ):
        raise ValueError("non-physical FV coefficients")
    face = 1.0 / (
        (edges[1:-1] - y[:-1]) / gamma[:-1] + (y[1:] - edges[1:-1]) / gamma[1:]
    )
    # Face conductance is Gamma_f / centre-centre distance and equals
    # 1/(d_owner/Gamma_owner + d_neighbor/Gamma_neighbor).
    wall = gamma[0] if wall_gamma is None else wall_gamma
    walltop = gamma[-1] if wall_gamma is None else wall_gamma
    lcoef = wall / (y[0] - edges[0])
    rcoef = walltop / (edges[-1] - y[-1])
    volumes = np.diff(edges)
    diagonal = np.zeros(n, float)
    diagonal[0] += lcoef
    diagonal[-1] += rcoef
    diagonal[:-1] += face
    diagonal[1:] += face
    diagonal += sink * volumes
    rhs = source * volumes
    rhs[0] += lcoef * bc0
    rhs[-1] += rcoef * bc1
    band = np.zeros((3, n))
    band[1] = diagonal
    band[0, 1:] = -face
    band[2, :-1] = -face
    sol = solve_banded(
        (1, 1), band, rhs, overwrite_ab=False, overwrite_b=False, check_finite=True
    )
    return sol, (diagonal, face, rhs)


def residual(f, operator):
    diag, off, rhs = operator
    r = diag * f - rhs
    r[:-1] -= off * f[1:]
    r[1:] -= off * f[:-1]
    return r


def state_terms(y, U, k, omega, c):
    d = np.minimum(y, c.height - y)
    Ugrad = np.gradient(U, y, edge_order=2)
    kg = np.gradient(k, y, edge_order=2)
    wg = np.gradient(omega, y, edge_order=2)
    cross_grad = kg * wg
    wi = np.maximum(omega, c.omega_floor)
    ki = np.maximum(k, c.k_floor)
    cd = np.maximum(2.0 * c.sigma_w2 / wi * cross_grad, 1e-10)
    arg1 = np.minimum(
        np.maximum(np.sqrt(ki) / (c.beta_star * wi * d), 500.0 * NU / (d * d * wi)),
        4.0 * c.sigma_w2 * ki / (cd * d * d),
    )
    arg2 = np.maximum(
        2.0 * np.sqrt(ki) / (c.beta_star * wi * d), 500.0 * NU / (d * d * wi)
    )
    f1 = np.tanh(np.minimum(arg1, 30.0) ** 4)
    f2 = np.tanh(np.minimum(arg2, 30.0) ** 2)
    nut = c.a1 * ki / np.maximum(c.a1 * wi, np.abs(Ugrad) * f2)
    production = np.minimum(
        nut * Ugrad * Ugrad, c.production_limiter * c.beta_star * ki * wi
    )
    sigk = f1 * c.sigma_k1 + (1.0 - f1) * c.sigma_k2
    sigw = f1 * c.sigma_w1 + (1.0 - f1) * c.sigma_w2
    beta = f1 * c.beta1 + (1.0 - f1) * c.beta2
    gamma = f1 * c.gamma1 + (1.0 - f1) * c.gamma2
    cross = 2.0 * (1.0 - f1) * c.sigma_w2 / wi * cross_grad
    return {
        "d": d,
        "grad": Ugrad,
        "f1": f1,
        "f2": f2,
        "nut": nut,
        "production": production,
        "sigk": sigk,
        "sigw": sigw,
        "beta": beta,
        "gamma": gamma,
        "cross": cross,
    }


_DEFAULT_PRECURSOR_CONTROLS = PrecursorControls()


def solve_precursor(c: PrecursorControls = _DEFAULT_PRECURSOR_CONTROLS):
    edges, y = grid(c)
    eta = y / c.height
    U = 6.0 * c.bulk_velocity * eta * (1.0 - eta)
    k = np.maximum(
        c.k_floor, 1.5 * (0.06 * c.bulk_velocity) ** 2 * np.sin(np.pi * eta) ** 2
    )
    omega = np.maximum(
        c.omega_floor,
        np.maximum(200.0, 6.0 * NU / (c.beta1 * np.minimum(y, c.height - y) ** 2)),
    )
    pressure_gradient = 0.0
    passed = False
    diagnostics = {}
    for it in range(1, c.max_iterations + 1):
        terms = state_terms(y, U, k, omega, c)
        omega_wall_l = 6.0 * NU / (c.beta1 * y[0] ** 2)
        omega_wall_r = 6.0 * NU / (c.beta1 * (c.height - y[-1]) ** 2)
        # Solve omega, then k, then U. Use local implicit destruction, not
        # explicit negative destruction that can violate the M-matrix.
        omega_prod = terms["gamma"] * np.where(
            terms["nut"] > 0, terms["production"] / np.maximum(terms["nut"], 1e-40), 0.0
        )
        source_w = omega_prod + np.maximum(terms["cross"], 0.0)
        sink_w = terms["beta"] * omega + np.maximum(-terms["cross"], 0.0) / np.maximum(
            omega, c.omega_floor
        )
        omega_solve, _opw = fv_diffusion(
            edges,
            y,
            NU + terms["sigw"] * terms["nut"],
            sink_w,
            source_w,
            omega_wall_l,
            omega_wall_r,
            wall_gamma=NU,
        )
        omega_solve = np.maximum(omega_solve, c.omega_floor)
        new_w = (1.0 - c.damping) * omega + c.damping * omega_solve
        # Updated omega and blending before reassembling k.
        updated = state_terms(y, U, k, new_w, c)
        k_solve, _opk = fv_diffusion(
            edges,
            y,
            NU + updated["sigk"] * updated["nut"],
            c.beta_star * new_w,
            updated["production"],
            0.0,
            0.0,
            wall_gamma=NU,
        )
        k_solve = np.maximum(k_solve, c.k_floor)
        new_k = (1.0 - c.damping) * k + c.damping * k_solve
        ku = state_terms(y, U, new_k, new_w, c)
        Ubase, _opU = fv_diffusion(
            edges,
            y,
            NU + ku["nut"],
            np.zeros_like(y),
            np.ones_like(y),
            0.0,
            0.0,
            wall_gamma=NU,
        )
        avg = np.sum(Ubase * np.diff(edges)) / c.height
        if not np.isfinite(avg) or avg <= 0.0:
            raise ArithmeticError("1D U solve has invalid bulk response")
        pressure_gradient = c.bulk_velocity / avg
        new_U = (1.0 - c.damping) * U + c.damping * pressure_gradient * Ubase
        change = max(
            float(np.max(np.abs(new_U - U)) / c.bulk_velocity),
            float(np.max(np.abs(new_k - k) / np.maximum(1.0, new_k))),
            float(np.max(np.abs(new_w - omega) / np.maximum(1.0, new_w))),
        )
        U, k, omega = new_U, new_k, new_w
        if not (
            np.isfinite(U).all()
            and np.isfinite(k).all()
            and np.isfinite(omega).all()
            and np.all(k > 0.0)
            and np.all(omega > 0.0)
        ):
            raise ArithmeticError("nonphysical 1D SST iterate")
        if it == 1 or it % 25 == 0 or change < 1e-6:
            # Nonlinear reassembled FV residual: frozen-source linear residual
            # alone would be invalid here. Include actual k*omega and cross.
            tt = state_terms(y, U, k, omega, c)
            _, oU = fv_diffusion(
                edges,
                y,
                NU + tt["nut"],
                np.zeros_like(y),
                np.full_like(y, pressure_gradient),
                0.0,
                0.0,
                wall_gamma=NU,
            )
            _, ok = fv_diffusion(
                edges,
                y,
                NU + tt["sigk"] * tt["nut"],
                c.beta_star * omega,
                tt["production"],
                0.0,
                0.0,
                wall_gamma=NU,
            )
            sw = (
                tt["gamma"]
                * np.where(
                    tt["nut"] > 0, tt["production"] / np.maximum(tt["nut"], 1e-40), 0.0
                )
                + tt["cross"]
            )
            _, ow = fv_diffusion(
                edges,
                y,
                NU + tt["sigw"] * tt["nut"],
                tt["beta"] * omega,
                sw,
                omega_wall_l,
                omega_wall_r,
                wall_gamma=NU,
            )

            def scaled_defect(f, op):
                r = residual(f, op)
                diag, off, rhs = op
                # Discrete rowwise residual normalized by maximum row action.
                scale = np.abs(diag * f) + np.abs(rhs)
                scale[:-1] += np.abs(off * f[1:])
                scale[1:] += np.abs(off * f[:-1])
                return float(np.max(np.abs(r) / np.maximum(scale, 1e-30)))

            diagnostics = {
                "iteration": it,
                "update_relative": change,
                "U_equation_defect": scaled_defect(U, oU),
                "k_equation_defect": scaled_defect(k, ok),
                "omega_equation_defect": scaled_defect(omega, ow),
                "bulk_velocity_m_s": float(np.dot(U, np.diff(edges)) / c.height),
                "mean_velocity_error": float(
                    abs(np.dot(U, np.diff(edges)) / c.height - c.bulk_velocity)
                    / c.bulk_velocity
                ),
                "momentum_pressure_source_m_s2": float(pressure_gradient),
                "minimum_k": float(k.min()),
                "minimum_omega": float(omega.min()),
                "maximum_nut_over_nu": float(tt["nut"].max() / NU),
                "omega_face_left": float(omega_wall_l),
                "omega_face_right": float(omega_wall_r),
            }
            if change <= c.tolerance and all(
                diagnostics[x] <= c.tolerance
                for x in (
                    "U_equation_defect",
                    "k_equation_defect",
                    "omega_equation_defect",
                    "mean_velocity_error",
                )
            ):
                passed = True
                break
    # Independent resolved-wall gate: smooth-wall low-Re SST should not
    # claim a verified precursor when the first cell is in the log layer.
    utau0 = math.sqrt(NU * max(U[0], 0.0) / y[0])
    utau1 = math.sqrt(NU * max(U[-1], 0.0) / (c.height - y[-1]))
    yplus0 = utau0 * y[0] / NU
    yplus1 = utau1 * (c.height - y[-1]) / NU
    diagnostics["first_cell_yplus_left"] = float(yplus0)
    diagnostics["first_cell_yplus_right"] = float(yplus1)
    diagnostics["cf_left"] = float(2.0 * utau0**2 / c.bulk_velocity**2)
    diagnostics["cf_right"] = float(2.0 * utau1**2 / c.bulk_velocity**2)
    diagnostics["discrete_equations_converged"] = bool(passed)
    diagnostics["wall_resolved"] = bool(max(yplus0, yplus1) <= 1.0)
    diagnostics["converged"] = bool(passed and diagnostics["wall_resolved"])
    diagnostics["configuration"] = asdict(c)
    return {
        "edges": edges,
        "y": y,
        "U": U,
        "k": k,
        "omega": omega,
        "diagnostics": diagnostics,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=192)
    parser.add_argument("--iterations", type=int, default=2500)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = solve_precursor(
        PrecursorControls(n_cells=args.cells, max_iterations=args.iterations)
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "y": result["y"].tolist(),
                "U": result["U"].tolist(),
                "k": result["k"].tolist(),
                "omega": result["omega"].tolist(),
                "diagnostics": result["diagnostics"],
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(result["diagnostics"], indent=2))
    if not result["diagnostics"]["converged"]:
        raise SystemExit("PRECURSOR_1D_NOT_CONVERGED: refuse physical profile export")


if __name__ == "__main__":
    main()
