#!/usr/bin/env python3
"""meshio -> CFDX-HDF5 mesh bridge."""
from __future__ import annotations
import argparse, collections, pathlib, sys
import h5py
import meshio
import numpy as np

VOLUME_FACES = {
 "tetra": ((0,2,1),(0,1,3),(1,3,2),(2,0,3)),
 "hexahedron": ((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)),
 "wedge": ((0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5)),
 "pyramid": ((0,1,2,3),(0,4,1),(1,2,4),(2,3,4),(3,0,4)),
 "voxel": ((0,4,6,2),(1,3,7,5),(0,1,5,4),(2,6,7,3),(0,2,3,1),(4,5,7,6)),
}
SURFACE = {"triangle":3,"quad":4,"triangle6":3,"quad8":4,"quad9":4}
CORNER_COUNT = {"tetra":4,"hexahedron":8,"wedge":6,"pyramid":5,"voxel":8,"tetra10":4,"hexahedron20":8,"hexahedron27":8,"wedge15":6,"wedge18":6,"pyramid13":5,"pyramid14":5,"triangle6":3,"quad8":4,"quad9":4}
BASE_TYPE = {"tetra10":"tetra","hexahedron20":"hexahedron","hexahedron27":"hexahedron","wedge15":"wedge","wedge18":"wedge","pyramid13":"pyramid","pyramid14":"pyramid"}

def patch_type(name):
    n=name.lower()
    if "inlet" in n:return 1
    if "outlet" in n:return 2
    if "symmetr" in n:return 3
    if "periodic" in n:return 4
    if "interface" in n:return 5
    if "empty" in n:return 6
    return 0 if "wall" in n else 7

def physical_names(mesh):
    out={}
    for name,value in getattr(mesh,"field_data",{}).items():
        try: out[int(value[0])]=name
        except (IndexError,TypeError,ValueError): pass
    return out

def block_tags(mesh,key,i):
    values=getattr(mesh,"cell_data",{}).get(key)
    return None if values is None or i>=len(values) else np.asarray(values[i])

def build(mesh):
    points=np.asarray(mesh.points,dtype=np.float64)
    if points.ndim!=2 or points.shape[1] not in (2,3): raise ValueError(f"unsupported point shape {points.shape}")
    if points.shape[1]==2: points=np.column_stack((points,np.zeros(len(points))))
    volumes=[]; surfaces=[]; skipped=collections.Counter()
    for bi,b in enumerate(mesh.cells):
        ctype=b.type; base=BASE_TYPE.get(ctype,ctype)
        if base in VOLUME_FACES:
            ncorner=CORNER_COUNT.get(ctype,len(VOLUME_FACES[base][0]))
            volumes += [(base,[int(x) for x in row[:ncorner]],bi) for row in np.asarray(b.data)]
        elif ctype in SURFACE:
            ncorner=SURFACE[ctype]
            surfaces += [(ctype,[int(x) for x in row[:ncorner]],bi) for row in np.asarray(b.data)]
        elif ctype=="polyhedron":
            for row in b.data: volumes.append(("polyhedron",[[int(x) for x in np.asarray(face)] for face in row],bi))
        else: skipped[ctype]+=len(b.data)
    if not volumes and not surfaces: raise ValueError("no supported cells")
    cells=volumes if volumes else surfaces
    face_map={}; faces=[]; owner=[]; neighbour=[]; cell_faces=[]
    for ci,(ctype,nodes,_) in enumerate(cells):
        if volumes and ctype=="polyhedron": templates=tuple(range(len(nodes)))
        elif volumes: templates=VOLUME_FACES[ctype]
        else: templates=tuple((i,(i+1)%len(nodes)) for i in range(len(nodes)))
        refs=[]
        for local in templates:
            face=nodes[local] if volumes and ctype=="polyhedron" else ([nodes[i] for i in local] if volumes else [nodes[local[0]],nodes[local[1]]])
            key=tuple(sorted(face)); fid=face_map.get(key)
            if fid is None:
                fid=len(faces); face_map[key]=fid; faces.append(face); owner.append(ci); neighbour.append(-1)
            elif neighbour[fid]!=-1: raise ValueError(f"non-manifold face {key}")
            else: neighbour[fid]=ci
            refs.append(fid)
        cell_faces.append(refs)
    return points,faces,owner,neighbour,cell_faces,skipped

def fnv1a_update(hash_value,array):
    for byte in np.asarray(array).tobytes(order="C"):
        hash_value ^= byte; hash_value=(hash_value*1099511628211)&0xFFFFFFFFFFFFFFFF
    return hash_value

def hash_hex(hash_value): return f"{hash_value:016x}"

def write(path,mesh,topo):
    points,faces,owner,neighbour,cell_faces,skipped=topo
    fv=np.asarray([v for f in faces for v in f],dtype=np.uint64)
    fo=np.zeros(len(faces)+1,dtype=np.uint64)
    for i,f in enumerate(faces): fo[i+1]=fo[i]+len(f)
    cf=np.asarray([f for c in cell_faces for f in c],dtype=np.uint64)
    co=np.zeros(len(cell_faces)+1,dtype=np.uint64)
    for i,c in enumerate(cell_faces): co[i+1]=co[i]+len(c)
    owner_a=np.asarray(owner,dtype=np.uint64); neighbour_a=np.asarray(neighbour,dtype=np.int64)
    lookup={tuple(sorted(f)):i for i,f in enumerate(faces)}
    names=physical_names(mesh); patches=collections.OrderedDict()
    for bi,b in enumerate(mesh.cells):
        if b.type not in SURFACE: continue
        tags=block_tags(mesh,"gmsh:physical",bi)
        for j,row in enumerate(np.asarray(b.data)):
            ncorner=SURFACE[b.type]; fid=lookup.get(tuple(sorted(int(x) for x in row[:ncorner])))
            if fid is None: continue
            pname=names.get(int(tags[j]),f"physical_{int(tags[j])}") if tags is not None and j<len(tags) else "boundary"
            patches.setdefault(pname,[]).append(fid)
    if not patches:
        ids=[i for i,n in enumerate(neighbour) if n<0]
        if ids: patches["boundary"]=ids
    # Every topological boundary face must belong to a patch. Physical
    # groups are optional/incomplete in real Gmsh files, so preserve their
    # named patches and explicitly collect any otherwise-unassigned boundary
    # faces in a generic patch instead of producing an invalid mesh.
    assigned={fid for values in patches.values() for fid in values}
    missing=[i for i,n in enumerate(neighbour) if n < 0 and i not in assigned]
    if missing:
        patches.setdefault("boundary", []).extend(missing)
        print("warning: unassigned boundary faces mapped to generic boundary patch: "
              + ",".join(str(i) for i in missing), file=sys.stderr)

    ids=[]; offsets=[0]; meta=[]
    for name,values in patches.items():
        values=sorted(set(values)); ids.extend(values); offsets.append(len(ids))
        start=values[0] if values else 0
        meta.append(f"{name}:{start}:{len(values)}:{patch_type(name)}")
    topology=1469598103934665603
    topology=fnv1a_update(topology,fv)
    if len(faces)>0: topology=fnv1a_update(topology,fo)
    topology=fnv1a_update(topology,owner_a); topology=fnv1a_update(topology,neighbour_a)
    topology=fnv1a_update(topology,cf)
    if len(cell_faces)>0: topology=fnv1a_update(topology,co)
    geometry=topology
    geometry=fnv1a_update(geometry,points[:,0].astype(np.float64,copy=False))
    geometry=fnv1a_update(geometry,points[:,1].astype(np.float64,copy=False))
    geometry=fnv1a_update(geometry,points[:,2].astype(np.float64,copy=False))
    with h5py.File(path,"w") as h:
        h.attrs["format"]="CFDX-HDF5-mesh-v1"
        h.attrs.create("format_version","1",dtype=h5py.string_dtype(encoding="ascii",length=2))
        h.attrs.create("schema_version","1",dtype=h5py.string_dtype(encoding="ascii",length=2))
        h.attrs.create("cfdx_version","0.7",dtype=h5py.string_dtype(encoding="ascii",length=4))
        h.attrs.create("topology_hash",hash_hex(topology),dtype=h5py.string_dtype(encoding="ascii",length=17))
        h.attrs.create("mesh_hash",hash_hex(geometry),dtype=h5py.string_dtype(encoding="ascii",length=17))
        h.attrs["n_points"]=str(len(points)); h.attrs["n_faces"]=str(len(faces)); h.attrs["n_cells"]=str(len(cell_faces))
        h.create_dataset("points",data=points); h.create_dataset("face_vertices",data=fv); h.create_dataset("face_offsets",data=fo)
        h.create_dataset("owner",data=owner_a); h.create_dataset("neighbour",data=neighbour_a)
        h.create_dataset("cell_faces",data=cf); h.create_dataset("cell_offsets",data=co)
        boundary_metadata=";".join(meta)
        if boundary_metadata: h.attrs.create("boundary_patches",boundary_metadata,dtype=h5py.string_dtype(encoding="ascii",length=len(boundary_metadata)+1))
        h.create_dataset("patch_face_ids",data=np.asarray(ids,dtype=np.uint64)); h.create_dataset("patch_face_offsets",data=np.asarray(offsets,dtype=np.uint64))
    if skipped: print("warning: skipped cell types: "+", ".join(f"{k}={v}" for k,v in sorted(skipped.items())),file=sys.stderr)

def main():
    p=argparse.ArgumentParser(); p.add_argument("input",type=pathlib.Path); p.add_argument("output",type=pathlib.Path); p.add_argument("--file-format")
    a=p.parse_args()
    try: mesh=meshio.read(a.input,file_format=a.file_format); write(a.output,mesh,build(mesh)); return 0
    except Exception as e: print(f"CFDX meshio import failed: {e}",file=sys.stderr); return 2
if __name__=="__main__": raise SystemExit(main())
