#!/usr/bin/env python3
"""Generate the reproducible VMFL036 Re=100 sphere mesh."""
import sys
import gmsh

if len(sys.argv) != 2:
    raise SystemExit("usage: generate_vmfl036_mesh.py OUTPUT.msh")

out = sys.argv[1]
gmsh.initialize()
try:
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("VMFL036_Re100")
    fluid = gmsh.model.occ.addCylinder(-10.0, 0.0, 0.0, 20.0, 0.0, 0.0, 10.0)
    sphere = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 0.5)
    cut, _ = gmsh.model.occ.cut([(3, fluid)], [(3, sphere)], removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()
    volumes = [tag for dim, tag in cut if dim == 3]
    if len(volumes) != 1:
        raise RuntimeError(f"expected one fluid volume, got {volumes}")

    groups = {"sphere": [], "inlet": [], "outlet": [], "farfield": []}
    for tag in (tag for dim, tag in gmsh.model.getEntities(2)):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        dx, dy, dz = xmax-xmin, ymax-ymin, zmax-zmin
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

    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.10)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 1.5)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.write(out)
    print("VMFL036_MESH: PASS nodes=" + str(len(gmsh.model.mesh.getNodes()[0])))
finally:
    gmsh.finalize()
