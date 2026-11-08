"""Bounded *actual polyhedral cell* tile reader for CFDX HDF5 volume meshes.

Unlike the boundary surface LOD, every cell in a requested contiguous range
is reconstructed from its original face connectivity. Only that range, the
faces it references, and its points are read. The caller may render one tile
at a time on CPU/GPU via VTK; no GPU-compute, global grid, streaming executive,
AMR or distributed rendering capability is implied by this reader.
"""
from __future__ import annotations

import math
import stat
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from .progressive_fields import _field


MAX_TILE_CELLS = 512
MAX_TILE_BYTES = 16 * 1024 * 1024


@dataclass(frozen=True, slots=True)
class VolumePartition:
    source: str
    source_identity: tuple[int, int, int, int]
    original_cell_ids: np.ndarray
    original_point_ids: np.ndarray
    points: np.ndarray
    polyhedron_stream: np.ndarray
    field_name: str | None
    field_values: np.ndarray | None
    total_cells: int
    max_bytes: int

    @property
    def fraction(self) -> float:
        return len(self.original_cell_ids) / self.total_cells

    def as_pyvista(self, pv):
        """VTK unstructured grid; GPU support depends on the active renderer."""
        types = np.full(len(self.original_cell_ids), int(pv.CellType.POLYHEDRON), dtype=np.uint8)
        grid = pv.UnstructuredGrid(self.polyhedron_stream, types, self.points)
        grid.cell_data["CFDX_original_cell_id"] = self.original_cell_ids
        grid.point_data["CFDX_original_point_id"] = self.original_point_ids
        if self.field_name is not None:
            grid.cell_data[self.field_name] = self.field_values
        return grid


def _dataset(h5: h5py.File, name: str, rank: int, kind: str) -> h5py.Dataset:
    link = h5.get(name, getlink=True)
    data = h5.get(name)
    if (not isinstance(link, h5py.HardLink) or not isinstance(data, h5py.Dataset)
            or data.is_virtual or data.external or data.ndim != rank
            or data.dtype.kind not in kind or data.dtype.itemsize > 8):
        raise ValueError(f"volume tile rejected unsafe HDF5 dataset: {name}")
    if data.chunks and math.prod(data.chunks) * data.dtype.itemsize > 8 * 1024 * 1024:
        raise ValueError(f"volume tile HDF5 compressed chunk exceeds budget: {name}")
    return data


def read_volume_partition(path: str|Path, start_cell: int, *, count: int = 128,
                          field: str|None = None, max_bytes: int = MAX_TILE_BYTES) -> VolumePartition:
    if (type(start_cell) is not int or start_cell < 0 or type(count) is not int
            or not 1 <= count <= MAX_TILE_CELLS
            or type(max_bytes) is not int or not 1024*1024 <= max_bytes <= MAX_TILE_BYTES):
        raise ValueError("invalid volume tile cell count/start or memory budget")
    path = Path(path)
    if path.is_symlink():
        raise ValueError("volume partition cannot read symbolic-link input")
    def identity():
        s = path.stat()
        if not stat.S_ISREG(s.st_mode):
            raise ValueError("volume partition source is not a regular file")
        return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns)
    before = identity()
    with h5py.File(path,"r") as h5:
        points = _dataset(h5,"points",2,"fiu")
        fv = _dataset(h5,"face_vertices",1,"iu")
        fo = _dataset(h5,"face_offsets",1,"iu")
        cf = _dataset(h5,"cell_faces",1,"iu")
        co = _dataset(h5,"cell_offsets",1,"iu")
        total = co.size - 1
        if (total < 1 or start_cell >= total or start_cell+count > total
                or points.shape[1:] != (3,) or points.shape[0] < 4):
            raise ValueError("volume tile range does not fit original mesh")
        # Reading offsets is O(count), NOT O(total mesh size).
        offsets = np.asarray(co[start_cell:start_cell+count+1],dtype=np.int64)
        if (offsets[0] < 0 or offsets[-1] > cf.size or np.any(np.diff(offsets) < 4)
                or int(offsets[-1]-offsets[0]) > count*64):
            raise ValueError("volume tile cell offsets/face count invalid")
        cell_faces = np.asarray(cf[int(offsets[0]):int(offsets[-1])],dtype=np.int64)
        unique_faces = np.unique(cell_faces)
        if (len(unique_faces) > count*64 or len(unique_faces)==0
                or unique_faces[0] < 0 or unique_faces[-1] >= fo.size-1):
            raise ValueError("volume tile references invalid faces")
        starts = np.asarray(fo[unique_faces],dtype=np.int64)
        stops = np.asarray(fo[unique_faces+1],dtype=np.int64)
        sizes = stops-starts
        if (np.any(sizes < 3) or np.any(sizes > 32) or np.any(starts < 0)
                or np.any(stops > fv.size) or int(sizes.sum()) > 65536):
            raise ValueError("volume tile has malformed or excessive face vertices")
        faces = {int(fid): np.asarray(fv[int(a):int(b)],dtype=np.int64)
                 for fid,a,b in zip(unique_faces,starts,stops)}
        all_ids = np.concatenate(list(faces.values()))
        unique_points = np.unique(all_ids)
        if len(unique_points) > 65536 or unique_points[0] < 0 or unique_points[-1] >= points.shape[0]:
            raise ValueError("volume tile has invalid original point IDs")
        coords = np.asarray(points[unique_points],dtype=np.float64)
        if not np.isfinite(coords).all():
            raise ValueError("volume tile has non-finite geometry")
        # np.searchsorted on sorted global IDs maps original face vertex IDs
        # to local VTK IDs without loading unrelated points.
        local_faces = {fid: np.searchsorted(unique_points,vertices)
                       for fid,vertices in faces.items()}
        stream: list[int] = []
        for local_cell in range(count):
            begin = int(offsets[local_cell]-offsets[0]); end = int(offsets[local_cell+1]-offsets[0])
            face_ids = cell_faces[begin:end]
            if len(set(map(int,face_ids))) != len(face_ids):
                raise ValueError("volume tile contains duplicate cell face IDs")
            # Avoid handing VTK a surface with missing/duplicated edges: a
            # polyhedron must have two incident faces per undirected edge.
            edges: dict[tuple[int,int],int] = {}
            vertices_in_cell: set[int] = set()
            stream.append(len(face_ids))
            for fid in face_ids:
                vertices = local_faces[int(fid)]
                if len(set(map(int,vertices))) != len(vertices):
                    raise ValueError("volume tile polyhedral face repeats vertices")
                for k in range(len(vertices)):
                    a,b = int(vertices[k]),int(vertices[(k+1)%len(vertices)])
                    key=(min(a,b),max(a,b))
                    edges[key]=edges.get(key,0)+1
                    vertices_in_cell.add(a)
                stream.append(len(vertices))
                stream.extend(map(int,vertices))
            if len(vertices_in_cell) < 4 or not edges or any(v != 2 for v in edges.values()):
                raise ValueError("volume tile has open/nonmanifold polyhedral cell")
        values = None
        if field is not None:
            ds = _field(h5,field)
            if ds.size != total:
                raise ValueError("volume tile scalar does not correspond to CFDX cell count")
            values = np.asarray(ds[start_cell:start_cell+count],dtype=np.float64)
            if not np.isfinite(values).all():
                raise ValueError("volume tile scalar contains non-finite data")
        cell_ids = np.arange(start_cell,start_cell+count,dtype=np.int64)
        stream_array = np.asarray(stream,dtype=np.int64)
        consumed = (offsets.nbytes+cell_faces.nbytes+starts.nbytes+stops.nbytes+all_ids.nbytes
                    +unique_points.nbytes+coords.nbytes+stream_array.nbytes
                    +sum(v.nbytes for v in faces.values()) + (0 if values is None else values.nbytes))
        if consumed > max_bytes:
            raise ValueError("volume tile exceeds memory budget")
    if before != identity():
        raise ValueError("volume tile source changed while reading")
    return VolumePartition(str(path.resolve()),before,cell_ids,unique_points,coords,
                           stream_array,field,values,total,max_bytes)


def main(argv: list[str]|None=None) -> int:
    import argparse
    import json
    p=argparse.ArgumentParser(description="Safely read CFDX polyhedral volume tiles")
    p.add_argument("mesh",type=Path)
    p.add_argument("--start",type=int,default=0)
    p.add_argument("--count",type=int,default=128)
    p.add_argument("--field")
    args=p.parse_args(argv)
    try:
        tile=read_volume_partition(args.mesh,args.start,count=args.count,field=args.field)
    except (ValueError,OSError) as exc:
        p.exit(2,f"Volume tile refused: {exc}\n")
    print(json.dumps({"mesh":tile.source,"total_cells":tile.total_cells,
                      "first_cell":int(tile.original_cell_ids[0]),
                      "cells_in_tile":len(tile.original_cell_ids),
                      "fraction":tile.fraction,"points_loaded":len(tile.points),
                      "renderer":"external_vtk_cpu_or_opengl_not_device_certified",
                      "scientific_convergence_certificate":False},indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
