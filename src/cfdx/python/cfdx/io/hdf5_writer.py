"""CFDX HDF5 case writer — writes case.cfdx.h5 from a ConversionResult.

This module bridges the Python conversion layer and the C++ thin-wrapper
adapters.  After the Python adapter populates a ``ConversionResult``
(mesh + setup + fields + GapAnalysis), this writer serialises it to the
CFDX HDF5 interchange format (``case.cfdx.h5``) that the C++ side reads
with ``read_case_cfdx_h5()``.

The mesh is stored in CFDX CSR (compressed-sparse-row) topology order so
that the C++ reader (``cfdx::io::read_mesh_hdf5``) can ingest it directly.
"""

from __future__ import annotations

import json
from typing import Optional

import h5py
import numpy as np

from cfdx.io.interfaces import ConversionResult
from cfdx.io.schema import CFDX_SCHEMA_VERSION


_ELEMENT_FACES: dict[str, list[list[int]]] = {
    "tetra": [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]],
    "hexahedron": [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
                   [1, 2, 6, 5], [2, 3, 7, 6], [0, 4, 7, 3]],
    "wedge": [[0, 2, 1], [3, 5, 4], [0, 1, 4, 3], [1, 2, 5, 4], [0, 3, 5, 2]],
    "pyramid": [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [0, 4, 3]],
    "triangle": [[0, 1, 2]],
    "quad": [[0, 1, 2, 3]],
}

_TYPE_ALIASES: dict[str, str] = {
    "tetrahedron": "tetra",
    "tetra": "tetra",
    "hexahedron": "hexahedron",
    "hex": "hexahedron",
    "wedge": "wedge",
    "prism": "wedge",
    "pyramid": "pyramid",
    "triangle": "triangle",
    "tri": "triangle",
    "quad": "quad",
    "quadrilateral": "quad",
}


def _meshio_to_cfdx(mesh: dict) -> Optional[dict]:
    """Convert a meshio-style mesh dict to CFDX CSR arrays.

    Returns ``None`` when the mesh contains no usable volumetric surface
    elements.
    """
    if not mesh or mesh.get("points") is None:
        return None

    points = np.asarray(mesh["points"], dtype=np.float64)
    if points.ndim == 1:
        points = points.reshape(-1, 3)
    n_points = points.shape[0]

    cells = mesh.get("cells", {})
    if not cells:
        return None

    cell_list: list[tuple[str, np.ndarray]] = []
    for cell_type, conn in cells.items():
        cfdx_type = _TYPE_ALIASES.get(cell_type, cell_type)
        if cfdx_type not in _ELEMENT_FACES:
            continue
        for i in range(len(conn)):
            cell_list.append((cfdx_type, np.asarray(conn[i], dtype=np.uint64)))

    if not cell_list:
        return None

    face_map: dict[tuple, list[tuple[int, int]]] = {}
    face_vertices: list[list[int]] = []
    cell_to_faces: list[list[int]] = [[] for _ in range(len(cell_list))]

    for cell_id, (cell_type, nodes) in enumerate(cell_list):
        for face_def in _ELEMENT_FACES[cell_type]:
            face = tuple(int(nodes[i]) for i in face_def)
            canonical = tuple(sorted(face))
            face_index = len(face_vertices)
            face_vertices.append(list(face))
            face_map.setdefault(canonical, []).append((cell_id, face_index))
            cell_to_faces[cell_id].append(face_index)

    n_faces = len(face_vertices)
    owner = np.zeros(n_faces, dtype=np.uint64)
    neighbour = np.full(n_faces, -1, dtype=np.int64)

    for occurrences in face_map.values():
        if len(occurrences) == 1:
            cell_id, face_index = occurrences[0]
            owner[face_index] = cell_id
            neighbour[face_index] = -1
        elif len(occurrences) >= 2:
            (c1, f1), (c2, f2) = occurrences[0], occurrences[1]
            owner[f1] = c1
            neighbour[f1] = c2
            owner[f2] = c2
            neighbour[f2] = c1
            for _, other_fi in occurrences[2:]:
                neighbour[other_fi] = -1

    fv_flat: list[int] = []
    for fv in face_vertices:
        fv_flat.extend(fv)

    cf_flat: list[int] = []
    for cf in cell_to_faces:
        cf_flat.extend(cf)

    face_offsets_arr = np.zeros(n_faces + 1, dtype=np.uint64)
    for i, fv in enumerate(face_vertices):
        face_offsets_arr[i + 1] = face_offsets_arr[i] + len(fv)

    cell_offsets_arr = np.zeros(len(cell_list) + 1, dtype=np.uint64)
    for i, cf in enumerate(cell_to_faces):
        cell_offsets_arr[i + 1] = cell_offsets_arr[i] + len(cf)

    return {
        "points": points.flatten(),
        "face_vertices": np.array(fv_flat, dtype=np.uint64),
        "face_offsets": face_offsets_arr,
        "owner": owner,
        "neighbour": neighbour,
        "cell_faces": np.array(cf_flat, dtype=np.uint64),
        "cell_offsets": cell_offsets_arr,
        "n_points": n_points,
        "n_faces": n_faces,
        "n_cells": len(cell_list),
    }


def write_case_cfdx_h5(result: ConversionResult, output_path) -> None:
    """Write a full ConversionResult to a ``case.cfdx.h5`` file.

    The file layout matches the C++ ``read_mesh_hdf5()`` reader:

    * Root-level datasets: ``points``, ``face_vertices``, ``face_offsets``,
      ``owner``, ``neighbour``, ``cell_faces``, ``cell_offsets``.
    * Root-level string attributes: ``source_solver``, ``source_format``,
      ``source_version``, ``source_case_path``, ``case_setup_json``,
      ``gap_report_json``.
    * Mesh-level attributes: ``boundary_patches`` (optional), integrity
      hashes.
    * ``/fields/scalar/<name>`` and ``/fields/vector/<name>`` datasets.
    """
    import pathlib
    output_path = pathlib.Path(output_path)

    with h5py.File(output_path, "w") as f:
        f.attrs["format_version"] = "1"
        f.attrs["schema_version"] = str(CFDX_SCHEMA_VERSION)
        f.attrs["cfdx_version"] = "0.7"

        source = result.source
        f.attrs["source_solver"] = source.solver or "unknown"
        f.attrs["source_format"] = source.format or "unknown"
        f.attrs["source_version"] = source.version or "unknown"
        f.attrs["source_case_path"] = source.case_path or ""
        f.attrs["source_case_name"] = source.case_name or ""

        if result.setup is not None:
            setup_json = result.setup.model_dump(mode="json")
            f.attrs["case_setup_json"] = json.dumps(setup_json, default=str)

        gap_json = _gap_report_to_json(result.gap_report)
        f.attrs["gap_report_json"] = gap_json

        if result.mesh is not None and result.mesh.get("points") is not None:
            cfdx_mesh = _meshio_to_cfdx(result.mesh)
            if cfdx_mesh is not None:
                _write_mesh_topology(f, cfdx_mesh)

        fields_grp = f.create_group("fields")
        if result.scalar_fields:
            sgrp = fields_grp.create_group("scalar")
            for name, data in result.scalar_fields:
                if data is not None and len(data) > 0:
                    sgrp.create_dataset(name, data=np.asarray(data, dtype=np.float64))
        if result.vec_fields:
            vgrp = fields_grp.create_group("vector")
            for name, data in result.vec_fields:
                if data is not None and len(data) > 0:
                    arr = np.asarray(data, dtype=np.float64)
                    if arr.ndim == 1:
                        arr = arr.reshape(-1, 1)
                    vgrp.create_dataset(name, data=arr)


def _write_mesh_topology(f: h5py.File, m: dict) -> None:
    f.create_dataset("points", data=m["points"])
    f.create_dataset("face_vertices", data=m["face_vertices"])
    f.create_dataset("face_offsets", data=m["face_offsets"])
    f.create_dataset("owner", data=m["owner"])
    f.create_dataset("neighbour", data=m["neighbour"])
    f.create_dataset("cell_faces", data=m["cell_faces"])
    f.create_dataset("cell_offsets", data=m["cell_offsets"])
    f.attrs["mesh_topology"] = "cfdx-csr-v1"


def _gap_report_to_json(gap) -> str:
    if gap is None:
        return json.dumps({
            "has_blocking": False,
            "summary": {
                "supported": 0, "approximated": 0,
                "unsupported_nonblocking": 0, "unsupported_blocking": 0,
                "unavailable": 0, "total": 0,
            },
            "findings": [],
        })

    findings = []
    for f in gap.findings:
        findings.append({
            "severity": f.severity.value,
            "category": f.category,
            "feature": f.feature,
            "detail": f.detail,
            "suggestion": f.suggestion,
        })

    return json.dumps({
        "has_blocking": gap.has_blocking(),
        "summary": {
            "supported": gap.n_supported(),
            "approximated": gap.n_approximated(),
            "unsupported_nonblocking": gap.n_unsupported_nonblocking(),
            "unsupported_blocking": gap.n_unsupported_blocking(),
            "unavailable": gap.n_unavailable(),
            "total": len(gap.findings),
        },
        "findings": findings,
    })
