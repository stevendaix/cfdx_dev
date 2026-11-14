#!/usr/bin/env python3
"""Experimental self-contained hybrid TH-03 initialisation on the TARGET mesh.

Builds a physically motivated wall-damped velocity and positive k/omega,
projects face mass flux through a sparse finite-volume potential solve, and
writes a NEW CFDX HDF5 case.  This is a preconditioned starting state, NOT
an exact Navier--Stokes, SST, or wall-function solution. Final solver gates
must not change. No coarse-grid solution or external VTU is required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

from scripts.generate_th03_openfoam_blockmesh import (
    K_IN,
    NU,
    OMEGA_IN,
    RHO,
    VELOCITY,
    H,
)


def _geometry(f: h5py.File) -> dict:
    xyz = np.asarray(f["points"][:], float)
    fo = np.asarray(f["face_offsets"][:], np.int64)
    fv = np.asarray(f["face_vertices"][:], np.int64)
    co = np.asarray(f["cell_offsets"][:], np.int64)
    cf = np.asarray(f["cell_faces"][:], np.int64)
    owner = np.asarray(f["owner"][:], np.int64)
    nb = np.asarray(f["neighbour"][:], np.int64)
    n = int(f.attrs["n_cells"])
    nf = len(owner)
    if len(co) != n + 1 or len(fo) != nf + 1 or len(nb) != nf:
        raise ValueError("inconsistent cell/face CSR geometry")
    centres = np.empty((n, 3))
    for i in range(n):
        face_ids = cf[co[i] : co[i + 1]]
        verts = np.unique(np.concatenate([fv[fo[j] : fo[j + 1]] for j in face_ids]))
        centres[i] = xyz[verts].mean(axis=0)
    face_centres = np.empty((nf, 3))
    sf = np.empty((nf, 3))
    for i in range(nf):
        verts = xyz[fv[fo[i] : fo[i + 1]]]
        if len(verts) < 3:
            raise ValueError("face with fewer than 3 vertices")
        face_centres[i] = verts.mean(axis=0)
        normal = 0.5 * np.sum(np.cross(verts, np.roll(verts, -1, axis=0)), axis=0)
        toward = (
            centres[nb[i]] - centres[owner[i]]
            if nb[i] >= 0
            else face_centres[i] - centres[owner[i]]
        )
        if np.dot(normal, toward) < 0:
            normal = -normal
        if np.linalg.norm(normal) <= 1e-19:
            raise ValueError("degenerate face normal")
        sf[i] = normal
    patch_ids = np.asarray(f["patch_face_ids"][:], np.int64)
    po = np.asarray(f["patch_face_offsets"][:], np.int64)
    meta = str(f.attrs["boundary_patches"]).split(";")
    if len(meta) != len(po) - 1:
        raise ValueError("boundary patch ordering mismatch")
    patches = {}
    for i, record in enumerate(meta):
        name, _start, count, _type = record.split(":")
        ids = patch_ids[po[i] : po[i + 1]]
        if len(ids) != int(count) or np.any(nb[ids] >= 0):
            raise ValueError("boundary patch topology mismatch")
        patches[name] = ids
    return {
        "centres": centres,
        "face_centres": face_centres,
        "sf": sf,
        "owner": owner,
        "nb": nb,
        "patches": patches,
    }


def _profile(
    xy: np.ndarray,
    wall_face_xy: np.ndarray,
    *,
    delta_h: float = 1.5,
    ramp_start_h: float = -110.0,
    ramp_end_h: float = -4.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bounded 1/7-power BL surrogate, with inlet uniform and mass per x-slice.

    This is *not* the Driver--Seegmiller measured inlet profile; tests must
    reject physical qualification without a separately verified profile.
    """
    if delta_h <= 0 or ramp_end_h <= ramp_start_h:
        raise ValueError("invalid BL ramp parameters")
    x, y = xy.T
    xh = x / H
    yh = y / H
    ylo = np.where(xh < 0.0, 1.0, 0.0)
    h = 9.0 - ylo
    if np.any((yh <= ylo) | (yh >= 9.0)):
        raise ValueError("cell centre outside physical channel")
    # Smooth step to avoid discontinuous acceleration at the startup patch.
    t = np.clip((xh - ramp_start_h) / (ramp_end_h - ramp_start_h), 0.0, 1.0)
    ramp = t * t * (3.0 - 2.0 * t)

    def transverse_factor(eta, low):
        d = np.minimum(eta - low, 9.0 - eta)
        return np.minimum(1.0, np.maximum(d, 0.0) / delta_h) ** (1.0 / 7.0)

    shape = (1.0 - ramp) + ramp * transverse_factor(yh, ylo)
    # Integrate the same 1/7-power profile to preserve inlet volumetric flow.
    eta = np.linspace(0.0, 1.0, 129)
    sample_y = ylo[:, None] + h[:, None] * eta[None, :]
    sample = (1.0 - ramp[:, None]) + ramp[:, None] * transverse_factor(
        sample_y, ylo[:, None]
    )
    integral = np.trapezoid(sample, eta, axis=1) * h
    ux = VELOCITY * 8.0 * shape / integral
    if not np.all(np.isfinite(ux)) or np.any(ux <= 0):
        raise ValueError("nonphysical initial velocity profile")
    # Actual nearest wall distance includes vertical step and both flat walls.
    wall_dist = cKDTree(wall_face_xy).query(xy, k=1)[0]
    if not np.all(np.isfinite(wall_dist)) or np.any(wall_dist <= 0):
        raise ValueError("invalid sampled wall distance")
    # Damped positive turbulence in near-wall cells, consistent with k -> 0.
    # omega is first-cell oriented and capped to prevent singular seeds.
    damper = np.minimum(1.0, wall_dist / (delta_h * H))
    k = np.maximum(1e-10, K_IN * damper**2)
    omega = np.maximum(OMEGA_IN, 6.0 * NU / (0.075 * wall_dist**2))
    omega = np.minimum(omega, 1.0e9)
    # Startup boundary is symmetry: do not invent a wall layer there.
    k[xh < ramp_start_h] = K_IN
    omega[xh < ramp_start_h] = OMEGA_IN
    return np.column_stack((ux, np.zeros(len(ux)), np.zeros(len(ux)))), k, omega


def _raw_flux(geo: dict, u: np.ndarray) -> np.ndarray:
    owner = geo["owner"]
    nb = geo["nb"]
    sf = geo["sf"]
    nf = len(owner)
    phi = np.zeros(nf)
    interior = nb >= 0
    phi[interior] = RHO * np.einsum(
        "ij,ij->i", 0.5 * (u[owner[interior]] + u[nb[interior]]), sf[interior]
    )
    # Boundary faces come from actual boundary contracts, not from ghosts.
    if "inlet" not in geo["patches"] or "outlet" not in geo["patches"]:
        raise ValueError("inlet and outlet patches required")
    inlet = geo["patches"]["inlet"]
    outlet = geo["patches"]["outlet"]
    phi[inlet] = RHO * VELOCITY * sf[inlet, 0]
    phi[outlet] = RHO * np.einsum("ij,ij->i", u[owner[outlet]], sf[outlet])
    # Everything else: wall, symmetry or front/back EMPTY => impermeable.
    return phi


def divergence(
    owner: np.ndarray, nb: np.ndarray, phi: np.ndarray, n: int
) -> np.ndarray:
    acc = np.bincount(owner, weights=phi, minlength=n)
    inside = nb >= 0
    acc -= np.bincount(nb[inside], weights=phi[inside], minlength=n)
    return acc


def project_mass_flux(geo: dict, raw: np.ndarray, n: int) -> tuple[np.ndarray, dict]:
    """Dirichlet-outlet FV potential projection preserving fixed boundary fluxes.

    Internal phi' = phi + c*(psi_owner - psi_neighbour); outlet psi_B=0.
    L psi = -D phi; L is SPD if every component is outlet-connected.
    """
    o = geo["owner"]
    nb = geo["nb"]
    sf = geo["sf"]
    cen = geo["centres"]
    fc = geo["face_centres"]
    inside = np.flatnonzero(nb >= 0)
    outlet = np.asarray(geo["patches"]["outlet"], dtype=int)
    norms = np.linalg.norm(sf, axis=1)
    dist = np.linalg.norm(cen[nb[inside]] - cen[o[inside]], axis=1)
    dist_out = np.linalg.norm(fc[outlet] - cen[o[outlet]], axis=1)
    if np.any(dist <= 0) or np.any(dist_out <= 0):
        raise ValueError("nonpositive projection distance")
    c_int = RHO * norms[inside] / dist
    c_out = RHO * norms[outlet] / dist_out
    ii = np.r_[o[inside], nb[inside], o[inside], nb[inside], o[outlet]]
    jj = np.r_[o[inside], nb[inside], nb[inside], o[inside], o[outlet]]
    data = np.r_[c_int, c_int, -c_int, -c_int, c_out]
    L = coo_matrix((data, (ii, jj)), shape=(n, n)).tocsc()
    rhs = -divergence(o, nb, raw, n)
    diag = L.diagonal()
    if np.any(diag <= 0) or not np.isfinite(diag).all():
        raise ValueError("nonpositive sparse projection diagonal")
    s = np.sqrt(diag)
    A = diags(1 / s) @ L @ diags(1 / s)
    factor = splu(A.tocsc())

    def solve(r):
        return factor.solve(r / s) / s

    psi = solve(rhs)
    # Iterative refinement on the original physical-scale linear equations.
    for _ in range(3):
        rr = rhs - L @ psi
        if np.max(np.abs(rr)) <= 1e-13 * max(1.0, np.max(np.abs(rhs))):
            break
        psi += solve(rr)
    phi = raw.copy()
    phi[inside] += c_int * (psi[o[inside]] - psi[nb[inside]])
    phi[outlet] += c_out * psi[o[outlet]]
    before = divergence(o, nb, raw, n)
    after = divergence(o, nb, phi, n)
    inlet = geo["patches"]["inlet"]
    btotal = max(1e-30, np.sum(np.abs(raw[inlet])))
    diag_report = {
        "before_l1_kg_s": float(np.sum(np.abs(before))),
        "after_l1_kg_s": float(np.sum(np.abs(after))),
        "normalized_after_l1": float(np.sum(np.abs(after)) / btotal),
        "after_linf_kg_s": float(np.max(np.abs(after))),
        "linear_backward_error": float(
            np.linalg.norm(L @ psi - rhs, np.inf)
            / max(np.linalg.norm(rhs, np.inf), 1e-30)
        ),
        "inlet_kg_s": float(np.sum(phi[inlet])),
        "outlet_kg_s": float(np.sum(phi[outlet])),
        "max_flux_change_kg_s": float(np.max(np.abs(phi - raw))),
        "boundary_immutable": bool(np.array_equal(phi[inlet], raw[inlet])),
        "n_cells": n,
    }
    # Gate on actual face divergence, NOT on the scalar solver's claim of convergence.
    if (
        diag_report["normalized_after_l1"] > 1.0e-7
        or not diag_report["boundary_immutable"]
    ):
        raise ArithmeticError(
            f"initial FV projection failed physical flux gate: {diag_report}"
        )
    return phi, diag_report


def initialise_case(
    source: Path,
    destination: Path,
    *,
    delta_h: float = 1.5,
    sst_style: str = "wall_aware",
    require_native_geometry: bool = True,
) -> dict:
    source = source.resolve(strict=True)
    if require_native_geometry:
        from scripts.require_th03_native_geometry import require_native_geometry as verify_native_geometry

        verify_native_geometry(source)
    destination = destination.resolve()
    if destination == source or destination.exists():
        raise ValueError("destination must be a new, absent HDF5 file")
    with h5py.File(source, "r") as f:
        if "fields" in f:
            raise ValueError("source case already has fields; refusing to modify")
        geo = _geometry(f)
    if sst_style not in (
        "wall_aware",
        "soft_wall",
        "uniform",
        "precursor_1d",
        "precursor_1d_u_only",
    ):
        raise ValueError("unsupported SST initialization style")
    n = len(geo["centres"])
    walls = np.concatenate([geo["patches"][p] for p in ("lowerWall", "upperWall")])
    u, k, omega = _profile(
        geo["centres"][:, :2], geo["face_centres"][walls, :2], delta_h=delta_h
    )
    precursor_evidence = None
    if sst_style in ("precursor_1d", "precursor_1d_u_only"):
        # This is a numerically solved, *fully developed channel* precursor,
        # not a claimed x-marching boundary layer or measured NASA inflow.
        # Keep the declared uniform inlet condition at x/H=-130; transition
        # the initial field only between -110 and the upstream step station.
        from scripts.solve_th03_precursor_1d import PrecursorControls, solve_precursor

        precursor = solve_precursor(PrecursorControls())
        precursor_evidence = precursor["diagnostics"]
        if not precursor_evidence["converged"]:
            raise ArithmeticError("unconverged 1D SST precursor cannot seed TH-03")
        yy = geo["centres"][:, 1] - H
        xx = geo["centres"][:, 0] / H
        upstream = xx < 0.0
        # First reach the fully-developed 1D channel state before the step.
        # Downstream, advect an initial memory of the incoming channel state;
        # a sharp reset to the old surrogate at x=0 introduces an artificial
        # k/omega jump unrelated to the actual separated-flow physics.
        t = np.clip((xx + 110.0) / 106.0, 0.0, 1.0)
        w = t * t * (3.0 - 2.0 * t)
        downstream = np.maximum(xx, 0.0)
        y_over_h = geo["centres"][:, 1] / H
        # The channel precursor is valid only for y/H in [1,9]; seed the
        # newly exposed lower cavity from the original post-step surrogate.
        bottom_gate = np.clip((y_over_h - 1.0) / 0.5, 0.0, 1.0)
        bottom_gate = bottom_gate * bottom_gate * (3.0 - 2.0 * bottom_gate)
        w[~upstream] = np.exp(-downstream[~upstream] / 8.0) * bottom_gate[~upstream]
        # Interpolate from the 1D FV cell-centre state. Reject extrapolated
        # stations, since they would silently break the 1D wall contract.
        if np.any(upstream & ((yy < 0.0) | (yy > 8.0 * H))):
            raise ValueError("1D precursor station lies outside upstream channel")
        sample = np.clip(yy, precursor["y"][0], precursor["y"][-1])
        ux = np.interp(sample, precursor["y"], precursor["U"])
        kval = np.interp(sample, precursor["y"], precursor["k"])
        wval = np.interp(sample, precursor["y"], precursor["omega"])
        u[:, 0] = (1.0 - w) * u[:, 0] + w * ux
        if sst_style == "precursor_1d":
            k[:] = (1.0 - w) * k + w * kval
            omega[:] = (1.0 - w) * omega + w * wval
        if not (
            np.isfinite(u).all()
            and np.isfinite(k).all()
            and np.isfinite(omega).all()
            and np.all(k > 0.0)
            and np.all(omega > 0.0)
        ):
            raise ArithmeticError("invalid interpolated precursor state")
    if sst_style == "soft_wall":
        # Reduce the severity of near-wall initial gradients, not the final SST BC.
        k = np.sqrt(k * K_IN)
        omega = np.sqrt(omega * OMEGA_IN)
    elif sst_style == "uniform":
        k = np.full(n, K_IN)
        omega = np.full(n, OMEGA_IN)
    phi_raw = _raw_flux(geo, u)
    phi, report = project_mass_flux(geo, phi_raw, n)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    try:
        with h5py.File(destination, "r+") as f:
            f.require_group("fields").require_group("vector").create_dataset(
                "U", data=u, compression="gzip"
            )
            s = f["fields"].require_group("scalar")
            s.create_dataset("k", data=k, compression="gzip")
            s.create_dataset("omega", data=omega, compression="gzip")
            f["fields"].require_group("face").create_dataset(
                "mdot", data=phi, compression="gzip"
            )
            # The velocity-only variant DOES NOT import certified 1D SST
            # fields: do not label its k/omega as a verified precursor.
            f.attrs["hybrid_initialisation"] = (
                "self_grid_channel_1d_U_k_omega_precursor_v2"
                if sst_style == "precursor_1d"
                else "self_grid_channel_1d_U_only_SST_surrogate_v2"
                if sst_style == "precursor_1d_u_only"
                else "self_grid_BL_SST_FV_flux_projection_v1"
            )
            f.attrs["hybrid_qualification"] = "NOT_QUALIFIED"
            f.attrs["hybrid_inlet_profile"] = (
                "CONVERGED_1D_CHANNEL_NOT_MEASURED_BFS_INFLOW"
                if sst_style in ("precursor_1d", "precursor_1d_u_only")
                else "UNVERIFIED_turbulent_BL_surrogate"
            )
    except BaseException:
        destination.unlink(missing_ok=True)
        raise
    report.update(
        {
            "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "output_case": str(destination),
            "fields_imported": ["U", "k", "omega", "mdot"],
            "velocity_range_m_s": [float(np.min(u[:, 0])), float(np.max(u[:, 0]))],
            "k_range_m2_s2": [float(np.min(k)), float(np.max(k))],
            "omega_range_s_inverse": [float(np.min(omega)), float(np.max(omega))],
            "initialisation": (
                "1D channel U/k/omega plus face-flux projection"
                if sst_style == "precursor_1d"
                else "1D channel U plus UNVERIFIED k/omega surrogate and face-flux projection"
                if sst_style == "precursor_1d_u_only"
                else "per-mesh BL surrogate plus sparse mass-flux projection"
            ),
            "sst_style": sst_style,
            "profile_validation": (
                "1D_U_K_OMEGA_DISCRETE_CERTIFIED_NOT_2D_QUALIFIED"
                if sst_style == "precursor_1d"
                else "1D_U_DISCRETE_CERTIFIED_K_OMEGA_UNVERIFIED_NOT_2D_QUALIFIED"
                if sst_style == "precursor_1d_u_only"
                else "NOT_QUALIFIED"
            ),
            "precursor_1d_evidence": precursor_evidence,
            "solver_flags": [
                "--initial-velocity-field",
                "--initial-sst-fields",
                "--initial-face-mass-flux",
                "--experimental-sst",
            ],
            "caveat": "Projected flux is a startup state, not a converged RANS solution; near-wall omega is approximate.",
        }
    )
    destination.with_suffix(".hybrid.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--allow-reconstructed-diagnostic",
        action="store_true",
        help="Explicitly allow reconstructed six-block geometry for diagnostic-only runs; never use for physical qualification.",
    )
    parser.add_argument("--delta-h", type=float, default=1.5)
    parser.add_argument(
        "--sst-style",
        choices=(
            "wall_aware",
            "soft_wall",
            "uniform",
            "precursor_1d",
            "precursor_1d_u_only",
        ),
        default="wall_aware",
    )
    args = parser.parse_args()
    print(
        json.dumps(
            initialise_case(
                args.case, args.output, delta_h=args.delta_h, sst_style=args.sst_style,
                require_native_geometry=not args.allow_reconstructed_diagnostic,
            ),
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
