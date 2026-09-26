#!/usr/bin/env python3
"""Generate the VMFL036 Re=100 3-D body-fitted sphere mesh."""

import math
import sys
import gmsh

if len(sys.argv) != 2:
    raise SystemExit("usage: generate_vmfl036_mesh.py OUTPUT.msh")

out = sys.argv[1]
gmsh.initialize()
try:
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("VMFL036_Re100")

    # Fluid cylinder: 20D long and 20D diameter, with the D=1 m sphere at x=0.
    fluid = gmsh.model.occ.addCylinder(-10.0, 0.0, 0.0, 20.0, 0.0, 0.0, 10.0)
    sphere = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 0.5)
    cut, _ = gmsh.model.occ.cut([(3, fluid)], [(3, sphere)], removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()

    volumes = [tag for dim, tag in cut if dim == 3]
    if len(volumes) != 1:
        raise RuntimeError(f"expected one fluid volume, got {volumes}")

    surfaces = [tag for dim, tag in gmsh.model.getEntities(2)]
    groups = {"sphere": [], "inlet": [], "outlet": [], "farfield": []}

    for tag in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
        if xmax < -9.999:
            groups["inlet"].append(tag)
        elif xmin > 9.999:
            groups["outlet"].append(tag)
        elif max(dx, dy, dz) <= 1.01:
            groups["sphere"].append(tag)
        else:
            groups["farfield"].append(tag)

    for name, tags in groups.items():
        if not tags:
            raise RuntimeError(f"no surfaces found for physical group {name}")
        pg = gmsh.model.addPhysicalGroup(2, tags)
        gmsh.model.setPhysicalName(2, pg, name)

    pv = gmsh.model.addPhysicalGroup(3, volumes)
    gmsh.model.setPhysicalName(3, pv, "fluid")

    # Resolve the sphere with ~10 elements across D and coarsen progressively away.
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.10)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 1.5)
    sphere_points = []
    for tag in groups["sphere"]:
        for dim, ptag in gmsh.model.getBoundary([(2, tag)], oriented=False, recursive=True):
            if dim == 0:
                sphere_points.append(ptag)
    for ptag in set(sphere_points):
        gmsh.model.mesh.setSize([(0, ptag)], 0.08)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.write(out)

    print(
        "VMFL036_MESH: PASS "
        f"nodes={gmsh.model.mesh.getNodes()[0].size} "
        f"surfaces={len(surfaces)} "
        f"sphere_faces={len(groups['sphere'])}"
    )
finally:
    gmsh.finalize()
