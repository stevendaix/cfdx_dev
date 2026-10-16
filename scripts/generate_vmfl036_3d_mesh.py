#!/usr/bin/env python3
"""Generate a reproducible 3-D spherical-shell mesh for the VMFL036 cross-check.

The sphere has diameter D=1 and the far-field sphere has radius 50 D. This
geometry intentionally mirrors the outer radius of the axisymmetric reference
case. The output is a Gmsh MSH 4.1 tetrahedral mesh with named physical groups:
fluid, sphere_wall and farfield.

Requires the Gmsh Python module. Mesh sizes are scaled by D and are deliberately
configurable so the coarse case can be used first to qualify the production
solver path before a refined campaign is attempted.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import gmsh


def generate(output: Path, lc_wall: float, lc_far: float, order: int = 1) -> None:
    if not (0.0 < lc_wall < lc_far):
        raise ValueError("Require 0 < lc_wall < lc_far")
    if order not in (1, 2):
        raise ValueError("Element order must be 1 or 2")

    diameter = 1.0
    sphere_radius = 0.5 * diameter
    farfield_radius = 50.0 * diameter
    output.parent.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.ElementOrder", order)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc_wall)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc_far)
        gmsh.model.add("vmfl036_sphere_3d")

        outer = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, farfield_radius)
        inner = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, sphere_radius)
        fluid, _ = gmsh.model.occ.cut([(3, outer)], [(3, inner)], removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()
        volumes = [tag for dim, tag in fluid if dim == 3]
        if len(volumes) != 1:
            raise RuntimeError(f"Expected one fluid volume, found {len(volumes)}")

        surface_entities = gmsh.model.getBoundary([(3, volumes[0])], oriented=False, recursive=False)
        wall_surfaces: list[int] = []
        farfield_surfaces: list[int] = []
        for dim, tag in surface_entities:
            if dim != 2:
                continue
            cx, cy, cz = gmsh.model.occ.getCenterOfMass(2, tag)
            radius = math.sqrt(cx * cx + cy * cy + cz * cz)
            if abs(radius - sphere_radius) < 1.0e-6 * diameter:
                wall_surfaces.append(tag)
            elif abs(radius - farfield_radius) < 1.0e-6 * diameter:
                farfield_surfaces.append(tag)

        if not wall_surfaces or not farfield_surfaces:
            raise RuntimeError(
                "Could not classify sphere wall and farfield surfaces "
                f"(wall={len(wall_surfaces)}, farfield={len(farfield_surfaces)})"
            )

        fluid_group = gmsh.model.addPhysicalGroup(3, volumes, 1)
        gmsh.model.setPhysicalName(3, fluid_group, "fluid")
        wall_group = gmsh.model.addPhysicalGroup(2, wall_surfaces, 2)
        gmsh.model.setPhysicalName(2, wall_group, "sphere_wall")
        farfield_group = gmsh.model.addPhysicalGroup(2, farfield_surfaces, 3)
        gmsh.model.setPhysicalName(2, farfield_group, "farfield")

        # Grade the mesh away from the sphere. The Threshold field is based on
        # distance to the sphere-wall surfaces, avoiding a uniform fine mesh
        # throughout the 50D far field.
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "FacesList", wall_surfaces)
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", lc_wall)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", lc_far)
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.05 * diameter)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", 10.0 * diameter)
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold)

        gmsh.model.mesh.generate(3)
        gmsh.write(str(output))
        print(
            "VMFL036_3D_MESH: PASS "
            f"output={output} lc_wall={lc_wall:g} lc_far={lc_far:g} "
            f"nodes={gmsh.model.mesh.getNodes()[0].size} "
            f"tetrahedra={sum(len(tags) for _, tags in gmsh.model.mesh.getElements(3)[1])}"
        )
    finally:
        gmsh.finalize()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lc-wall", type=float, default=0.10,
                        help="target element size at sphere wall, in units of D")
    parser.add_argument("--lc-far", type=float, default=4.0,
                        help="target element size in the far field, in units of D")
    parser.add_argument("--order", type=int, choices=(1, 2), default=1)
    args = parser.parse_args()
    generate(args.output, args.lc_wall, args.lc_far, args.order)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
