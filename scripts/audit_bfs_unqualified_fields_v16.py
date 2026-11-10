"""Audit CFDX BFS VTK cell fields without promoting an unconverged state.

This is a diagnostic of dimensional extrema and turbulence magnification, not
an experimental validation or a nonlinear residual certificate. Only ASCII
VTU CellData is handled; no interpolation or fabricated wall shear is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

REQUIRED = ("k", "omega", "nu_t", "u_x", "u_y", "u_z", "p", "wall_distance")


def audit_vtu(path: Path, *, u_ref: float, molecular_nu: float,
              inlet_k: float, inlet_omega: float) -> dict:
    for name, value in {"u_ref": u_ref, "molecular_nu": molecular_nu,
                        "inlet_k": inlet_k, "inlet_omega": inlet_omega}.items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    root = ET.parse(path).getroot()
    if root.attrib.get("type") != "UnstructuredGrid":
        raise ValueError("expected VTK UnstructuredGrid")
    pieces = root.findall("./UnstructuredGrid/Piece")
    if len(pieces) != 1:
        raise ValueError("exactly one VTU piece is required")
    piece = pieces[0]
    cells = int(piece.attrib["NumberOfCells"])
    if cells <= 0:
        raise ValueError("empty cell field")
    container = piece.find("CellData")
    if container is None:
        raise ValueError("missing CellData")
    fields = {}
    for element in container.findall("DataArray"):
        name = element.attrib.get("Name")
        if name not in REQUIRED:
            continue
        if name in fields or element.attrib.get("format") != "ascii" or element.attrib.get("NumberOfComponents", "1") != "1":
            raise ValueError("duplicate, non-scalar or non-ASCII field")
        values = [float(token) for token in (element.text or "").split()]
        if len(values) != cells or any(not math.isfinite(value) for value in values):
            raise ValueError("non-finite or mismatched cell field")
        fields[name] = values
    missing = set(REQUIRED) - set(fields)
    if missing:
        raise ValueError("missing CFDX fields: " + ", ".join(sorted(missing)))
    if any(fields["k"][i] < 0 or fields["nu_t"][i] < 0 or
           fields["omega"][i] <= 0 or fields["wall_distance"][i] <= 0
           for i in range(cells)):
        raise ValueError("nonphysical turbulence field or wall distance")
    speed = [math.sqrt(fields["u_x"][i] ** 2 + fields["u_y"][i] ** 2 +
                       fields["u_z"][i] ** 2) for i in range(cells)]
    if any(not math.isfinite(v) for v in speed):
        raise ValueError("overflowed velocity norm")
    def range_summary(name: str, values: list[float]) -> dict:
        high = max(range(len(values)), key=values.__getitem__)
        low = min(range(len(values)), key=values.__getitem__)
        return {"min": values[low], "min_cell": low,
                "max": values[high], "max_cell": high}
    ranges = {name: range_summary(name, values) for name, values in fields.items()}
    ranges["speed"] = range_summary("speed", speed)
    kmax = ranges["k"]["max"]
    result = {
        "schema": "CFDX-BFS-CELL-FIELDS/v16",
        "snapshot_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "cells": cells,
        "reference": {"u_ref_m_s": u_ref, "nu_m2_s": molecular_nu,
                      "inlet_k_m2_s2": inlet_k, "inlet_omega_s_inv": inlet_omega},
        "fields": ranges,
        "diagnostic_ratios": {
            "max_k_to_inlet_k": kmax / inlet_k,
            "max_omega_to_inlet_omega": ranges["omega"]["max"] / inlet_omega,
            "max_nu_t_to_nu": ranges["nu_t"]["max"] / molecular_nu,
            "max_speed_to_u_ref": ranges["speed"]["max"] / u_ref,
            "max_k_based_intensity_normalized_by_u_ref":
                math.sqrt(2.0 * kmax / 3.0) / u_ref,
            "negative_u_x_cells": sum(value < 0 for value in fields["u_x"]),
        },
        "algebraic_convergence": "NOT_CERTIFIED_BY_VTU",
        "physical_validation": "NOT_QUALIFIED",
        "reason": "requires committed U,p,k,omega conservation gates, native wall functions, 3-grid GCI and NASA Cf/Cp/profiles",
    }
    if any(not math.isfinite(v) for v in result["diagnostic_ratios"].values()):
        raise ValueError("dimensionless diagnostic overflow")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtu", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--u-ref", type=float, required=True)
    parser.add_argument("--nu", type=float, required=True)
    parser.add_argument("--inlet-k", type=float, required=True)
    parser.add_argument("--inlet-omega", type=float, required=True)
    args = parser.parse_args()
    result = audit_vtu(args.vtu, u_ref=args.u_ref, molecular_nu=args.nu,
                       inlet_k=args.inlet_k, inlet_omega=args.inlet_omega)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
