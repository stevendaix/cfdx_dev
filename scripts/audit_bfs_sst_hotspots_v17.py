"""Independent BFS SST hotspot audit on an accepted VTU cell snapshot.

No inference of nonlinear convergence or experimental agreement is made here.
The SST production limiter is evaluated from the returned cell k, omega,
nu_t and strain-rate; this is a local diagnostic, not a FV budget certificate.
Static pressure is p* - 2*rho*k/3 for CFDX SST pressure conventions.
"""
from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

SCHEMA = "CFDX-BFS-SST-HOTSPOTS/v17"
FIELDS = ("k", "omega", "nu_t", "sst_strain_rate", "p", "u_x", "u_y", "u_z")


def _values(tag: ET.Element) -> list[float]:
    if tag.get("format", "ascii") != "ascii" or tag.text is None:
        raise ValueError("Only uncompressed ASCII VTU arrays can be audited")
    values = [float(s) for s in tag.text.split()]
    if not values or not all(math.isfinite(v) for v in values):
        raise ValueError("VTU contains nonfinite/empty array")
    return values


def read_vtu(path: Path) -> tuple[dict[str, list[float]], list[list[float]]]:
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError("VTU UnstructuredGrid/Piece missing")
    n = int(piece.attrib["NumberOfCells"])
    arrays = piece.findall("./CellData/DataArray")
    fields = {item.attrib.get("Name", ""): _values(item) for item in arrays}
    if any(name not in fields or len(fields[name]) != n for name in FIELDS):
        raise ValueError("VTU is missing a required SST cell field")
    points_tag = piece.find("./Points/DataArray")
    if points_tag is None:
        raise ValueError("VTU points missing")
    points_flat = _values(points_tag)
    if len(points_flat) % 3:
        raise ValueError("VTU point coordinates malformed")
    points = [points_flat[j:j+3] for j in range(0, len(points_flat), 3)]
    cell_arrays = {item.attrib.get("Name", ""): item for item in piece.findall("./Cells/DataArray")}
    if "connectivity" not in cell_arrays or "offsets" not in cell_arrays:
        raise ValueError("VTU cell topology missing")
    indices = _values(cell_arrays["connectivity"])
    offsets = _values(cell_arrays["offsets"])
    if len(offsets) != n:
        raise ValueError("VTU cell offsets malformed")
    centers: list[list[float]] = []
    prev = 0
    for offset_float in offsets:
        end = int(offset_float)
        if end <= prev or end != offset_float or end > len(indices):
            raise ValueError("VTU invalid offsets")
        ids = indices[prev:end]
        if any(int(i) != i or i < 0 or i >= len(points) for i in ids):
            raise ValueError("VTU invalid point index")
        centers.append([math.fsum(points[int(i)][d] for i in ids)/len(ids) for d in range(3)])
        prev = end
    if prev != len(indices):
        raise ValueError("VTU connectivity not fully accounted for")
    return fields, centers


def audit(path: Path, rho: float, nu: float, u_ref: float,
          k_inlet: float, production_limiter: float = 10.0,
          beta_star: float = 0.09, top_n: int = 8) -> dict:
    if any(not math.isfinite(v) or v <= 0 for v in
           (rho, nu, u_ref, k_inlet, production_limiter, beta_star)) or top_n < 1:
        raise ValueError("invalid physical reference parameters")
    fields, centers = read_vtu(path)
    n = len(fields["k"])
    pstar = fields["p"]
    k = fields["k"]
    omega = fields["omega"]
    nut = fields["nu_t"]
    strain = fields["sst_strain_rate"]
    if any(not math.isfinite(v) or v < 0 for name in ("k", "nu_t", "sst_strain_rate")
           for v in fields[name]) or any(not math.isfinite(v) or v <= 0 for v in omega):
        raise ValueError("invalid turbulence field in snapshot")
    q_ref = 0.5 * rho * u_ref * u_ref

    def sample(i: int) -> dict:
        Pk_unlimited = rho * nut[i] * strain[i] ** 2
        Pk_limit = production_limiter * beta_star * rho * k[i] * omega[i]
        Pk = min(Pk_unlimited, Pk_limit)
        Dk = rho * beta_star * k[i] * omega[i]
        p_static = pstar[i] - 2 * rho * k[i] / 3
        result = {
            "cell": i, "centroid_m": centers[i], "k_m2_s2": k[i],
            "omega_s_inv": omega[i], "nut_m2_s": nut[i], "strain_s_inv": strain[i],
            "Pk_unlimited_W_m3": Pk_unlimited, "Pk_limited_W_m3": Pk,
            "k_destruction_W_m3": Dk, "net_k_source_W_m3": Pk-Dk,
            "production_limit_active": Pk_unlimited > Pk_limit,
            "production_to_destruction": Pk / Dk if Dk > 0 else None,
            "nut_to_nu": nut[i] / nu,
            "k_to_inlet": k[i] / k_inlet,
            "k_based_intensity_over_inlet_speed": math.sqrt(2 * k[i] / 3) / u_ref,
            "p_star_Pa": pstar[i], "p_static_Pa": p_static,
            "Cp_star_minus_Cp_static": (pstar[i] - p_static) / q_ref,
            "speed_m_s": math.sqrt(sum(fields[d][i] ** 2 for d in ("u_x", "u_y", "u_z"))),
        }
        if any(isinstance(v, float) and not math.isfinite(v) for v in result.values()):
            raise ValueError("overflow in SST hotspot diagnostic")
        return result

    tracked: dict[int, list[str]] = {}
    for name in ("k", "nu_t", "omega", "sst_strain_rate"):
        for i in sorted(range(n), key=lambda i: fields[name][i], reverse=True)[:top_n]:
            tracked.setdefault(i, []).append(name)
    cells = [dict(sample(i), ranked_in=tracked[i]) for i in sorted(tracked)]
    return {
        "schema": SCHEMA, "source": path.name, "cells": n,
        "status": "DIAGNOSTIC_ONLY_NOT_QUALIFIED",
        "pressure_convention": "p_static=p_star-(2/3)*rho*k; Cp requires experimental reference pressure offset",
        "reference": {"density_kg_m3": rho, "nu_m2_s": nu,
                      "inlet_speed_m_s": u_ref, "inlet_k_m2_s2": k_inlet,
                      "sst_production_limiter": production_limiter, "beta_star": beta_star},
        "hotspots": cells,
        "notice": "Production and destruction are local SST estimates; no face flux balances or convergence proven.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtu", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rho", type=float, default=1.0)
    parser.add_argument("--nu", type=float, default=1.56e-5)
    parser.add_argument("--u-ref", type=float, default=44.2)
    parser.add_argument("--k-inlet", type=float, default=0.00109)
    parser.add_argument("--production-limiter", type=float, default=10.0)
    parser.add_argument("--beta-star", type=float, default=0.09)
    parser.add_argument("--top", type=int, default=8)
    args = parser.parse_args()
    summary = audit(args.vtu, args.rho, args.nu, args.u_ref, args.k_inlet,
                    args.production_limiter, args.beta_star, args.top)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "cells": summary["cells"],
                     "hotspots": len(summary["hotspots"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
