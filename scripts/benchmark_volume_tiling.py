"""Bounded synthetic HDF5 volume tile benchmark (I/O only, not GPU FPS).

The benchmark uses independently defined hexahedral cells, preserves original
cell identity and reports wall time/RSS. It is not a physical CFD qualification
nor a benchmark of distributed VTK, GPU rendering or a 50M-cell production mesh.
"""
from __future__ import annotations

import argparse
import json
import resource
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np

# Execute with: PYTHONPATH=src/cfdx/python python -m scripts.benchmark_volume_tiling
from cfdx.volume_partitions import MAX_TILE_BYTES, read_volume_partition


_FACES=np.asarray(((0,3,2,1),(4,5,6,7),(0,1,5,4),
                   (3,7,6,2),(0,4,7,3),(1,2,6,5)),dtype=np.int64)


def write_synthetic_grid(path: Path, n: int):
    if not 1 <= n <= 100_000:
        raise ValueError("synthetic benchmark size outside 1..100000")
    x=np.arange(n,dtype=np.float64)
    xyz=np.zeros((n,8,3),dtype=np.float64)
    xyz[:,:,0]=x[:,None]+np.asarray([0,1,1,0,0,1,1,0])
    xyz[:,:,1]=np.asarray([0,0,1,1,0,0,1,1])
    xyz[:,:,2]=np.asarray([0,0,0,0,1,1,1,1])
    vertices=(_FACES[None,:,:]+(8*np.arange(n))[:,None,None]).reshape(-1)
    with h5py.File(path,"w") as h5:
        h5.create_dataset("points",data=xyz.reshape(-1,3),chunks=(min(8192,n*8),3))
        h5.create_dataset("face_vertices",data=vertices,chunks=(min(16384,n*24),))
        h5.create_dataset("face_offsets",data=np.arange(0,n*24+1,4,dtype=np.int64))
        h5.create_dataset("cell_faces",data=np.arange(n*6,dtype=np.int64))
        h5.create_dataset("cell_offsets",data=np.arange(0,n*6+1,6,dtype=np.int64))
        h5.create_dataset("fields/scalar/T",data=np.linspace(273.,373.,n))
        h5["fields/scalar/T"].attrs["units"]="K"


def run_sizes(sizes: list[int], tile: int):
    data=[]
    with tempfile.TemporaryDirectory(prefix="cfdx-tiles-") as root:
        for n in sizes:
            path=Path(root)/f"mesh_{n}.h5"
            write_synthetic_grid(path,n)
            start=max(0,n//2-tile//2)
            start=min(start,n-tile)
            begin=time.perf_counter()
            part=read_volume_partition(path,start,count=tile,field="/fields/scalar/T")
            seconds=time.perf_counter()-begin
            resident=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            if len(part.original_cell_ids)!=tile or part.original_cell_ids[0]!=start:
                raise RuntimeError("tile reader failed original IDs")
            data.append({"synthetic_cells":n,"tile_cells":tile,"seconds":round(seconds,6),
                         "maxrss_kib":resident,"tile_fractions":part.fraction,
                         "points_in_tile":len(part.points),"max_tile_bytes":MAX_TILE_BYTES,
                         "gpu_render_fps":None,"physical_case":False})
    return data


def main(argv=None):
    parser=argparse.ArgumentParser()
    parser.add_argument("--sizes",default="1000,10000")
    parser.add_argument("--tile",type=int,default=128)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args(argv)
    sizes=[int(x) for x in args.sizes.split(",")]
    if not sizes or any(n < args.tile or n > 100_000 for n in sizes):
        parser.error("all synthetic sizes must satisfy tile <= size <= 100000")
    data={"kind":"CFDX-VOLUME-TILE-IO-BENCHMARK","scope":"synthetic CPU I/O only",
          "complete_gpu_certification":False,"results":run_sizes(sizes,args.tile)}
    payload=json.dumps(data,indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload,encoding="utf-8")
    else:
        print(payload,end="")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
