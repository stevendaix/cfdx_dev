"""Inspect prescribed-velocity farfield face fluxes on a matched CFDX HDF5 mesh.

This is a diagnostic, NEVER a convergence/force/GCI certificate. A globally
zero signed mass defect does not mean that each prescribed-velocity boundary
face satisfies its kinematic flux condition.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import h5py

REQUIRED = ("face", "owner", "neighbour", "cx", "cy", "cz",
            "sx", "sy", "sz", "mass_flux")


def _patch_faces(mesh_file: Path, patch: str) -> tuple[set[int], object, object, dict]:
    with h5py.File(mesh_file, "r") as mesh:
        metadata = str(mesh.attrs["boundary_patches"])
        names = [s.split(":", 1)[0] for s in metadata.split(";")]
        if len(set(names)) != len(names) or patch not in names:
            raise ValueError("unknown or duplicated boundary patch")
        i = names.index(patch)
        offsets = mesh["patch_face_offsets"][:]
        if len(offsets) != len(names) + 1:
            raise ValueError("mesh patch offset count mismatch")
        ids = mesh["patch_face_ids"][int(offsets[i]):int(offsets[i + 1])]
        face_ids = [int(f) for f in ids]
        if len(face_ids) != len(set(face_ids)) or not face_ids:
            raise ValueError("empty or duplicate patch face ids")
        owners, neighs = mesh["owner"][:], mesh["neighbour"][:]
        # Boundary neighbour IDs are -1; tolerate an unsigned legacy encoding.
        identity_ncells = int(mesh.attrs["n_cells"])
        for f in face_ids:
            if f < 0 or f >= len(owners):
                raise ValueError("invalid patch face id")
            if 0 <= int(neighs[f]) < identity_ncells:
                raise ValueError("a patch contains an internal face")
        identity = {
            "n_cells": int(mesh.attrs["n_cells"]),
            "n_faces": int(mesh.attrs["n_faces"]),
            "mesh_hash": str(mesh.attrs.get("mesh_hash", "")),
            "topology_hash": str(mesh.attrs.get("topology_hash", "")),
        }
        if identity["n_faces"] != len(owners):
            raise ValueError("mesh face count mismatch")
        return set(face_ids), owners, neighs, identity


def analyze(mesh_file: Path, flux_csv: Path, *, patch: str = "farfield",
            velocity: tuple[float, float, float] = (1.0, 0.0, 0.0),
            density: float = 1.0) -> tuple[dict, list[dict]]:
    if not math.isfinite(density) or density <= 0 or any(
            not math.isfinite(u) for u in velocity):
        raise ValueError("density and prescribed velocity must be finite")
    selected, owners, neighs, identity = _patch_faces(mesh_file, patch)
    seen = set()
    samples = []
    face_count = 0
    sum_error = 0.0
    abs_error = 0.0
    abs_reference = 0.0
    reference_area = 0.0
    largest = (-1.0, -1)
    with flux_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("missing/duplicate CSV columns")
        if not set(REQUIRED).issubset(reader.fieldnames):
            raise ValueError("missing face-flux CSV columns")
        for lineno, row in enumerate(reader, start=2):
            if None in row or any(v is None or not v for v in row.values()):
                raise ValueError(f"invalid face-flux CSV row {lineno}")
            fid_text = row["face"]
            if not fid_text.isdecimal():
                raise ValueError(f"noninteger face id on line {lineno}")
            fid = int(fid_text)
            face_count += 1
            if fid in seen or fid != face_count - 1 or fid >= len(owners):
                raise ValueError(f"duplicated, missing, or reordered face id {fid}")
            seen.add(fid)
            owner_text = row["owner"]
            neighbour_text = row["neighbour"]
            if not owner_text.isdecimal():
                raise ValueError(f"invalid owner on line {lineno}")
            if not (neighbour_text.isdecimal() or (neighbour_text.startswith("-") and
                    neighbour_text[1:].isdecimal())):
                raise ValueError(f"invalid neighbour on line {lineno}")
            if int(owner_text) != int(owners[fid]):
                raise ValueError(f"face owner mismatch at {fid}")
            neighbour = int(neighbour_text)
            actual_neighbour = int(neighs[fid])
            if neighbour != (-1 if actual_neighbour < 0 or actual_neighbour >= identity["n_cells"] else actual_neighbour):
                raise ValueError(f"face neighbour mismatch at {fid}")
            if fid not in selected:
                continue
            xyz_area_flux = [float(row[k]) for k in ("cx", "cy", "cz", "sx", "sy", "sz", "mass_flux")]
            if not all(math.isfinite(v) for v in xyz_area_flux):
                raise ValueError(f"nonfinite patch face {fid}")
            x, y, z, sx, sy, sz, mass_flux = xyz_area_flux
            area = math.sqrt(sx*sx + sy*sy + sz*sz)
            if area <= 0:
                raise ValueError(f"degenerate area on patch face {fid}")
            expected = density * (velocity[0]*sx + velocity[1]*sy + velocity[2]*sz)
            delta = mass_flux - expected
            relative_to_face_speed = abs(delta) / (density * area * math.sqrt(sum(v*v for v in velocity))) if any(velocity) else 0.0
            reference_area += area
            sum_error += delta
            abs_error += abs(delta)
            abs_reference += abs(expected)
            samples.append({"face": fid, "owner": int(owner_text), "cx": x, "cy": y, "cz": z,
                            "area": area, "mass_flux": mass_flux, "prescribed_flux": expected,
                            "mismatch": delta, "mismatch_rel_face": relative_to_face_speed})
            if abs(delta) > largest[0]:
                largest = (abs(delta), fid)
    if face_count != identity["n_faces"] or len(selected) != len(samples):
        raise ValueError("face CSV incomplete or selected patch not fully represented")
    worst = next(s for s in samples if s["face"] == largest[1])
    result = {
        "status": "DIAGNOSTIC_ONLY_NOT_PHYSICAL_CERTIFICATE",
        "mesh": identity, "patch": patch, "selected_face_count": len(samples),
        "density": density, "prescribed_velocity": list(velocity),
        "mismatch_sum_signed": sum_error, "mismatch_l1": abs_error,
        "prescribed_flux_l1": abs_reference,
        "mismatch_l1_relative_to_prescribed_flux_l1": abs_error / max(abs_reference, 1e-300),
        "mismatch_linf": largest[0], "mismatch_worst_face": largest[1],
        "mismatch_worst_relative_to_face_speed": worst["mismatch_rel_face"],
        "worst_face": worst,
        "compatible_to_absolute_1e-10": largest[0] <= 1e-10,
        "note": "Signed sum can cancel even when individual faces violate prescribed velocity mass flux",
    }
    return result, sorted(samples, key=lambda s: abs(s["mismatch"]), reverse=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--flux", required=True, type=Path)
    parser.add_argument("--patch", default="farfield")
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--ux", type=float, default=1.0)
    parser.add_argument("--uy", type=float, default=0.0)
    parser.add_argument("--uz", type=float, default=0.0)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worst-output", type=Path)
    args = parser.parse_args()
    result, face_rows = analyze(args.mesh, args.flux, patch=args.patch,
                                density=args.density, velocity=(args.ux, args.uy, args.uz))
    content = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8")
    else:
        print(content)
    if args.worst_output:
        with args.worst_output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(face_rows[0]))
            writer.writeheader()
            writer.writerows(face_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
