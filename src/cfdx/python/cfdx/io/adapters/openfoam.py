"""OpenFOAM case directory adapter for CFDX I/O.

Parses OpenFOAM case directory structure:
  - constant/polyMesh/points, faces, owner, neighbour, boundary
  - 0/ field files (U, p, T, etc.)

Corresponds to the C++ header:
  src/cfdx/io/openfoam/openfoam_importer.h
"""

from __future__ import annotations

import json
import re
from pathlib import Path


import numpy as np

from cfdx.io.interfaces import ConversionResult, SolverAdapter
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.schema import BCType, BCValueType, BoundarySpec, CaseSetup, SourceInfo


_FIELD_DIMS = {"scalar": 1, "vector": 3, "symmTensor": 6, "tensor": 9}


def _extract_foam_list(content: str):
    """Extract the first parenthesised list block from a FoamFile ascii file.

    Returns ``(lines, count)`` where ``lines`` are the comment-stripped,
    non-empty lines inside the block and ``count`` is the integer written just
    before the opening parenthesis (``None`` when absent).  Returns
    ``(None, None)`` when no list block is present.
    """
    content = re.sub(r"/\*.*?\*/", " ", content, flags=re.DOTALL)
    content = re.sub(r"//[^\n]*", " ", content)
    start = content.find("(")
    if start < 0:
        return None, None

    depth = 0
    end = -1
    for i in range(start, len(content)):
        ch = content[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                end = i
                break
    if end < 0:
        return None, None

    count = None
    prefix_match = re.search(r"(\d+)\s*$", content[:start])
    if prefix_match:
        count = int(prefix_match.group(1))

    lines = []
    for raw in content[start + 1:end].splitlines():
        line = re.sub(r"/\*.*?\*/", " ", raw)
        line = re.sub(r"//.*$", "", line).strip()
        if line:
            lines.append(line)
    return lines, count


class OpenFOAMAdapter(SolverAdapter):
    solver_name = "OpenFOAM"
    format_name = "openfoam"

    def __init__(self, rules: MappingRules | None = None) -> None:
        self.rules = rules or MappingRules.load_default()
        self.setup = CaseSetup(source=SourceInfo(solver="OpenFOAM", format="openfoam"))

    def convert(self, case_path: str | Path, result: ConversionResult) -> bool:
        path = Path(str(case_path))
        
        if not path.is_dir():
            result.gap_report.unsupported_blocking(
                "case_path",
                "openfoam_case_directory",
                f"OpenFOAM case must be a directory: {path}",
            )
            return False

        poly_mesh = path / "constant" / "polyMesh"
        if not poly_mesh.exists():
            result.gap_report.unsupported_blocking(
                "case_path",
                "openfoam_polymesh",
                f"Missing constant/polyMesh directory: {poly_mesh}",
            )
            return False

        # Parse mesh
        self._mesh_parse_error = "unknown parser error"
        points, cells = self._parse_mesh(poly_mesh)
        if points is None or cells is None:
            result.gap_report.unsupported_blocking(
                "mesh", "openfoam_polymesh_parse",
                "Failed to parse OpenFOAM mesh: "
                f"{getattr(self, '_mesh_parse_error', 'unknown parser error')}"
            )
            return False

        result.mesh = {"points": points, **cells}
        cfdx_csr = cells.get("cfdx_csr", {})
        if cfdx_csr:
            self.setup.mesh_info.n_vertices = int(cfdx_csr["n_points"])
            self.setup.mesh_info.n_faces = int(cfdx_csr["n_faces"])
            self.setup.mesh_info.n_cells = int(cfdx_csr["n_cells"])
            self.setup.mesh_info.n_boundary_faces = int(
                np.count_nonzero(np.asarray(cfdx_csr["neighbour"]) < 0)
            )
            self.setup.mesh_info.n_internal_faces = (
                self.setup.mesh_info.n_faces - self.setup.mesh_info.n_boundary_faces
            )
            self.setup.mesh_info.n_patches = len(
                cfdx_csr.get("boundary_patch_assignment", [])
            )
            self.setup.mesh_info.has_polyhedral = True
            self.setup.mesh_info.cell_types = ["polyhedron"]
        result.gap_report.supported("mesh", "native OpenFOAM polyhedral CSR topology")

        # Parse boundary conditions
        self._parse_boundary(poly_mesh, result)
        result.gap_report.supported("boundary", "boundary patches")

        # Parse fields from 0/ directory
        zero_dir = path / "0"
        if zero_dir.exists():
            n_parsed = self._parse_fields(zero_dir, result)
            if n_parsed:
                result.gap_report.supported(
                    "fields", "0/ directory fields", f"{n_parsed} field(s) parsed"
                )
            else:
                result.gap_report.unavailable(
                    "fields",
                    "zero_directory",
                    f"No readable fields in {zero_dir}",
                )
        else:
            result.gap_report.unavailable(
                "fields", "zero_directory", "0/ directory not found"
            )

        self._translate_setup_dictionaries(path, result)
        if any(name == "T" for name, _ in result.scalar_fields):
            self.setup.energy_model = "energy"

        result.setup = self.setup
        result.source.case_path = str(path)
        result.source.solver = "OpenFOAM"
        result.source.format = "openfoam"

        return result.success

    def _translate_setup_dictionaries(
        self, case_path: Path, result: ConversionResult
    ) -> None:
        """Preserve setup dictionaries and map only explicitly supported controls."""
        from cfdx.io.openfoam_case_reader import read_openfoam_dictionaries

        try:
            dictionaries = read_openfoam_dictionaries(case_path)
        except (OSError, ValueError) as exc:
            # Keep source configuration recoverable when a dictionary uses
            # syntax outside the native parser's supported ASCII subset.
            raw_dictionaries = {}
            for relative in (
                "system/controlDict",
                "system/fvSchemes",
                "system/fvSolution",
                "system/fvOptions",
                "constant/transportProperties",
                "constant/physicalProperties",
                "constant/thermophysicalProperties",
                "constant/turbulenceProperties",
                "constant/radiationProperties",
                "constant/g",
            ):
                source_file = case_path / relative
                if source_file.is_file():
                    try:
                        raw_dictionaries[relative] = source_file.read_text()
                    except OSError:
                        continue
            self.setup.source_metadata["openfoam_raw_dictionaries_json"] = json.dumps(
                raw_dictionaries, sort_keys=True
            )
            result.gap_report.unsupported_nonblocking(
                "configuration",
                "openfoam_dictionary_parser",
                f"Setup dictionaries were preserved as raw text but not parsed: {exc}",
                "The native parser supports common ASCII dictionaries; inspect preserved raw text for unsupported syntax.",
            )
            return

        self.setup.source_metadata["openfoam_dictionaries_json"] = json.dumps(
            dictionaries, sort_keys=True, default=str
        )
        control = dictionaries.get("system/controlDict", {})
        application = str(control.get("application", "")).strip()
        if application:
            self.setup.source_metadata["openfoam_application"] = application

        def numeric_setting(mapping: dict, key: str) -> float | None:
            value = mapping.get(key)
            if isinstance(value, (int, float)):
                return float(value)
            try:
                return float(str(value))
            except (TypeError, ValueError):
                return None

        delta_t = numeric_setting(control, "deltaT")
        end_time = numeric_setting(control, "endTime")
        if delta_t is not None and delta_t > 0:
            self.setup.time_step = delta_t
        if end_time is not None:
            self.setup.end_time = end_time

        fv_schemes = dictionaries.get("system/fvSchemes", {})
        ddt_schemes = fv_schemes.get("ddtSchemes", {})
        ddt_default = str(ddt_schemes.get("default", "")).strip().lower()
        if ddt_default:
            self.setup.transient = ddt_default != "steadystate"
            self.setup.solver_mode = "transient" if self.setup.transient else "steady"
        elif application.lower() in {
            "pimplefoam", "icofoam", "interfoam", "rhopimplefoam",
            "rhointerfoam", "reactingfoam", "buoyantpimplefoam",
        }:
            self.setup.transient = True
            self.setup.solver_mode = "transient"

        if self.setup.transient and self.setup.time_step > 0 and self.setup.end_time > 0:
            self.setup.max_time_steps = int(self.setup.end_time / self.setup.time_step)

        turbulence = dictionaries.get("constant/turbulenceProperties", {})
        simulation_type = str(turbulence.get("simulationType", "")).strip()
        if simulation_type.lower() == "laminar":
            self.setup.turbulence_model = "laminar"
        elif simulation_type in {"RAS", "LES"}:
            model_dict = turbulence.get(simulation_type, {})
            model_name = model_dict.get("model", model_dict.get("RASModel", model_dict.get("LESModel", "")))
            if model_name:
                self.setup.turbulence_model = str(model_name)

        radiation = dictionaries.get("constant/radiationProperties", {})
        radiation_model = radiation.get("radiationModel")
        if radiation_model is not None:
            self.setup.radiation_model = str(radiation_model)

        for key in ("system/fvSchemes", "system/fvSolution"):
            if key in dictionaries:
                self.setup.numerics.raw_settings[f"openfoam.{key}"] = json.dumps(
                    dictionaries[key], sort_keys=True, default=str
                )

        if any(key in dictionaries for key in ("system/fvSchemes", "system/fvSolution")):
            result.gap_report.unsupported_nonblocking(
                "numerics",
                "openfoam_numerics_not_fully_mapped",
                "fvSchemes and fvSolution were preserved, but are not yet translated "
                "to equivalent CFDX numerical registry selections.",
            )
        result.gap_report.supported(
            "configuration",
            "openfoam_dictionary_snapshot",
            f"Preserved {len(dictionaries)} parsed OpenFOAM setup dictionaries.",
        )

    def import_results(self, results_path: str, result: ConversionResult) -> bool:
        """Import cell-centred scalar/vector fields from an OpenFOAM time.

        The latest numeric time directory is selected. The native parser supports
        ASCII OpenFOAM fields and reports binary or unsupported syntax explicitly.
        """
        path = Path(results_path).expanduser()
        if not path.is_dir():
            result.gap_report.unsupported_blocking(
                "results", "openfoam_results_path",
                f"OpenFOAM results path must be a case directory: {path}",
            )
            return False

        try:
            from cfdx.io.openfoam_case_reader import (
                read_openfoam_time_fields,
            )

            owner = self._parse_label_list(path / "constant" / "polyMesh" / "owner")
            cell_count = max(owner) + 1 if owner else None
            fields = read_openfoam_time_fields(path, "latest", cell_count=cell_count)
        except (
            FileNotFoundError,
            NotADirectoryError,
            ValueError,
        ) as exc:
            result.gap_report.unsupported_blocking(
                "results", "openfoam_results_read", str(exc),
            )
            return False

        if not fields:
            result.gap_report.unavailable(
                "results", "openfoam_empty_time",
                f"No readable fields found in the latest time directory of {path}",
            )
            return False

        imported = 0
        # Latest-time data supersedes any initial values parsed from 0/.
        result.scalar_fields = []
        result.vec_fields = []
        latest_time = next(iter(fields.values()))["time"]
        time_key = f"{latest_time:g}"
        result.available_fields[time_key] = sorted(fields)
        result.field_data_source = f"OpenFOAM time {time_key}"

        for name, field_data in fields.items():
            source_class = field_data["class"]
            values = np.asarray(field_data["internal"], dtype=np.float64)

            # Surface fields and point fields have different locations and
            # cannot be placed in CFDX's cell-centred field arrays.
            if not source_class.startswith("vol"):
                result.gap_report.unsupported_nonblocking(
                    "results", name,
                    f"Field class {source_class} is not a volume field; values were not imported.",
                )
                continue

            if source_class.endswith("ScalarField"):
                if values.ndim != 1:
                    result.gap_report.unsupported_nonblocking(
                        "results", name,
                        f"Expected a cell-centred scalar array; got shape {values.shape}.",
                    )
                    continue
                result.scalar_fields.append((name, values))
                imported += 1
            elif source_class.endswith("VectorField"):
                if values.ndim == 1 and values.size == 3:
                    values = values.reshape(1, 3)
                if values.ndim != 2 or values.shape[1] != 3:
                    result.gap_report.unsupported_nonblocking(
                        "results", name,
                        f"Expected a cell-centred 3-component vector array; got shape {values.shape}.",
                    )
                    continue
                result.vec_fields.append((name, values))
                imported += 1
            else:
                result.gap_report.unsupported_nonblocking(
                    "results", name,
                    f"Field class {source_class} is not yet representable in ConversionResult.",
                    "Tensor and other field types remain available in the source OpenFOAM case.",
                )

        if imported:
            result.gap_report.supported(
                "results", "openfoam_time_fields",
                f"Imported {imported} scalar/vector volume field(s) from time {time_key}.",
            )
        else:
            result.gap_report.unavailable(
                "results", "openfoam_no_supported_fields",
                f"No scalar/vector volume fields could be imported from time {time_key}.",
            )
        return imported > 0

    def _parse_mesh(self, poly_mesh: Path) -> tuple[np.ndarray | None, dict | None]:
        """Parse OpenFOAM polyMesh files."""
        try:
            for list_name in ("points", "faces", "owner", "neighbour"):
                list_path = poly_mesh / list_name
                if list_path.is_file():
                    header = list_path.read_bytes()[:4096].decode("ascii", errors="ignore")
                    if re.search(r"\bformat\s+binary\s*;", header):
                        raise ValueError(
                            f"Binary OpenFOAM polyMesh file is not supported: {list_path}"
                        )
            # Parse points
            points = self._parse_points(poly_mesh / "points")
            if points is None:
                return None, None

            # Parse faces
            faces = self._parse_faces(poly_mesh / "faces")
            if faces is None:
                return None, None

            # Parse owner
            owner = self._parse_label_list(poly_mesh / "owner")
            if owner is None:
                return None, None

            # Parse neighbour
            neighbour = self._parse_label_list(poly_mesh / "neighbour")

            # Parse boundary
            patches = self._parse_boundary_file(poly_mesh / "boundary")

            # Build cell connectivity from owner/neighbour
            cells = self._build_cells_from_faces(faces, owner, neighbour, patches)

            return points, cells

        except (OSError, ValueError, TypeError, IndexError, KeyError, UnicodeError) as e:
            self._mesh_parse_error = str(e)
            return None, None

    def _parse_points(self, points_file: Path) -> np.ndarray | None:
        """Parse an ASCII OpenFOAM points list."""
        try:
            content = points_file.read_text(encoding="utf-8")
            if re.search(r"\bformat\s+binary\s*;", content):
                raise ValueError("Binary OpenFOAM polyMesh is not supported by the native Python reader")

            data_lines, _count = _extract_foam_list(content)
            if not data_lines:
                return None

            points = []
            for line in data_lines:
                match = re.match(
                    r"^\s*(?:vector\s*)?\(?\s*([-\d.eE+,\s]+?)\s*\)?\s*;?\s*$", line
                )
                if not match:
                    continue
                coords = match.group(1).replace(",", " ").split()
                if len(coords) >= 3:
                    points.append([float(c) for c in coords[:3]])

            if not points:
                return None
            return np.array(points, dtype=np.float64)
        except (OSError, ValueError, TypeError, IndexError, KeyError, UnicodeError):
            return None

    def _parse_faces(self, faces_file: Path) -> list | None:
        """Parse OpenFOAM faces file.

        Handles both the compact ``4(0 1 2 3)`` form written by OpenFOAM and
        the parenthesised ``(0 1 2 3)`` form.
        """
        try:
            content = faces_file.read_text(encoding="utf-8")
            if re.search(r"\bformat\s+binary\s*;", content):
                raise ValueError("Binary OpenFOAM polyMesh is not supported by the native Python reader")

            data_lines, _ = _extract_foam_list(content)
            if data_lines is None:
                return None

            faces = []
            for line in data_lines:
                match = re.match(r"^\s*(\d*)\s*\(?\s*([-\d.eE+,\s]+?)\s*\)?\s*;?\s*$", line)
                if not match:
                    continue
                n_declared = int(match.group(1)) if match.group(1) else None
                indices = [int(tok) for tok in match.group(2).replace(",", " ").split()]
                if n_declared is not None and n_declared != len(indices):
                    continue
                if len(indices) < 3:
                    continue
                faces.append(indices)

            if not faces:
                return None
            return faces
        except (OSError, ValueError, TypeError, IndexError, KeyError, UnicodeError):
            return None

    def _parse_label_list(self, label_file: Path) -> list | None:
        """Parse OpenFOAM label list (owner, neighbour)."""
        try:
            content = label_file.read_text(encoding="utf-8")
            if re.search(r"\bformat\s+binary\s*;", content):
                raise ValueError("Binary OpenFOAM polyMesh is not supported by the native Python reader")

            data_lines, count = _extract_foam_list(content)
            if data_lines is None:
                return None

            labels = []
            for line in data_lines:
                token = line.strip().rstrip(';')
                if not token or token.startswith(("/*", "//")):
                    continue
                try:
                    labels.append(int(token))
                except ValueError:
                    return None

            # A leading count token is only dropped when the number of parsed
            # entries does not match it (i.e. the count was still in the list).
            if labels and count is not None and labels[0] == count == len(labels) - 1:
                labels = labels[1:]

            if not labels:
                return []
            return labels
        except (OSError, ValueError, TypeError, IndexError, KeyError, UnicodeError):
            return None

    def _parse_boundary_file(self, boundary_file: Path) -> list:
        """Parse named OpenFOAM boundary patches and their face ranges."""
        try:
            content = boundary_file.read_text()
        except OSError:
            return []

        content = re.sub(r"/\*.*?\*/", " ", content, flags=re.DOTALL)
        content = re.sub(r"//[^\n]*", " ", content)
        patch_pattern = re.compile(
            r"([A-Za-z_][A-Za-z0-9_.-]*)\s*\{([^}]*)\}", re.DOTALL
        )
        patches = []
        for match in patch_pattern.finditer(content):
            name, body = match.groups()
            type_match = re.search(r"\btype\s+([^;]+);", body)
            faces_match = re.search(r"\bnFaces\s+(\d+)\s*;", body)
            start_match = re.search(r"\bstartFace\s+(\d+)\s*;", body)
            if type_match is None or faces_match is None or start_match is None:
                continue
            patches.append(
                {
                    "name": name,
                    "type": type_match.group(1).strip(),
                    "n_faces": int(faces_match.group(1)),
                    "start_face": int(start_match.group(1)),
                }
            )
        return patches

    def _build_cells_from_faces(self, faces, owner, neighbour, patches):
        """Build native CFDX face/cell CSR topology from OpenFOAM polyMesh.

        OpenFOAM already stores the exact polyhedral topology as faces, owner,
        neighbour, and boundary patch ranges. Preserve that topology directly;
        do not pretend boundary surface triangles are volume cells.
        """
        if len(owner) != len(faces):
            raise ValueError(
                f"owner has {len(owner)} entries but faces has {len(faces)}"
            )
        if neighbour is None:
            raise ValueError("Cannot parse OpenFOAM neighbour list")
        if len(neighbour) > len(faces):
            raise ValueError(
                f"neighbour has {len(neighbour)} entries but faces has {len(faces)}"
            )

        n_faces = len(faces)
        neighbour_full = np.full(n_faces, -1, dtype=np.int64)
        if neighbour:
            neighbour_full[: len(neighbour)] = np.asarray(neighbour, dtype=np.int64)

        max_cell = max(
            max(owner, default=-1),
            int(np.max(neighbour_full)) if np.any(neighbour_full >= 0) else -1,
        )
        n_cells = max_cell + 1
        if n_cells <= 0:
            raise ValueError("OpenFOAM mesh contains no cells")

        cell_faces: list[list[int]] = [[] for _ in range(n_cells)]
        for face_id, owner_cell in enumerate(owner):
            if owner_cell < 0 or owner_cell >= n_cells:
                raise ValueError(f"Invalid owner cell {owner_cell} for face {face_id}")
            cell_faces[owner_cell].append(face_id)
            neighbour_cell = int(neighbour_full[face_id])
            if neighbour_cell >= 0:
                if neighbour_cell >= n_cells:
                    raise ValueError(
                        f"Invalid neighbour cell {neighbour_cell} for face {face_id}"
                    )
                cell_faces[neighbour_cell].append(face_id)

        face_offsets = np.zeros(n_faces + 1, dtype=np.uint64)
        face_vertex_ids: list[int] = []
        for face_id, face in enumerate(faces):
            if len(face) < 3:
                raise ValueError(
                    f"OpenFOAM face {face_id} has fewer than three vertices"
                )
            if any(vertex < 0 or vertex >= len(self._mesh_points) for vertex in face):
                raise ValueError(f"Invalid vertex index in OpenFOAM face {face_id}")
            face_vertex_ids.extend(int(vertex) for vertex in face)
            face_offsets[face_id + 1] = len(face_vertex_ids)

        cell_offsets = np.zeros(n_cells + 1, dtype=np.uint64)
        cell_face_ids: list[int] = []
        for cell_id, cell_face_list in enumerate(cell_faces):
            if not cell_face_list:
                raise ValueError(f"OpenFOAM cell {cell_id} has no faces")
            cell_face_ids.extend(cell_face_list)
            cell_offsets[cell_id + 1] = len(cell_face_ids)

        patch_assignment = []
        for patch in patches:
            first = int(patch["start_face"])
            count = int(patch["n_faces"])
            stop = first + count
            if first < len(neighbour) or first < 0 or stop > n_faces:
                raise ValueError(
                    f"Invalid boundary range for patch {patch['name']!r}: "
                    f"startFace={first}, nFaces={count}"
                )
            patch_assignment.append(
                {
                    "name": patch["name"],
                    "type": patch["type"],
                    "face_ids": list(range(first, stop)),
                }
            )

        boundary_face_ids = set(np.flatnonzero(neighbour_full < 0).tolist())
        assigned_face_ids = {
            face_id
            for patch in patch_assignment
            for face_id in patch["face_ids"]
        }
        if assigned_face_ids != boundary_face_ids:
            missing = sorted(boundary_face_ids - assigned_face_ids)
            extra = sorted(assigned_face_ids - boundary_face_ids)
            raise ValueError(
                "OpenFOAM boundary patch ranges do not cover exactly the boundary "
                f"faces (missing={missing[:8]}, extra={extra[:8]})"
            )

        cfdx_csr = {
            "points": np.asarray(self._mesh_points, dtype=np.float64),
            "face_vertices": np.asarray(face_vertex_ids, dtype=np.uint64),
            "face_offsets": face_offsets,
            "owner": np.asarray(owner, dtype=np.uint64),
            "neighbour": neighbour_full,
            "cell_faces": np.asarray(cell_face_ids, dtype=np.uint64),
            "cell_offsets": cell_offsets,
            "n_points": len(self._mesh_points),
            "n_faces": n_faces,
            "n_cells": n_cells,
            "boundary_patch_assignment": patch_assignment,
        }

        # Keep a surface-only representation for lightweight VTK export. The
        # HDF5 writer must use cfdx_csr above for actual solver topology.
        triangles = []
        for face in faces:
            for index in range(1, len(face) - 1):
                triangles.append([face[0], face[index], face[index + 1]])
        surface_cells = {
            "triangle": np.asarray(triangles, dtype=np.int64)
        } if triangles else {}

        return {"cells": surface_cells, "cfdx_csr": cfdx_csr}

    def _parse_boundary(self, poly_mesh: Path, result: ConversionResult):
        """Parse boundary patches and add to setup."""
        patches = self._parse_boundary_file(poly_mesh / "boundary")
        for patch in patches:
            bc_type = self._map_patch_type(patch['type'])
            bc_value = BCValueType.FIXED if bc_type in (BCType.INLET, BCType.OUTLET) else BCValueType.ZERO_GRADIENT
            
            self.setup.boundary_conditions.append(BoundarySpec(
                patch_name=patch['name'],
                type=bc_type,
                value_type=bc_value,
            ))

    def _map_patch_type(self, of_type: str) -> BCType:
        """Map OpenFOAM patch type to CFDX BCType."""
        mapping = {
            'patch': BCType.WALL,
            'wall': BCType.WALL,
            'inlet': BCType.INLET,
            'outlet': BCType.OUTLET,
            'pressureInletOutletVelocity': BCType.PRESSURE_OUTLET,
            'pressureOutlet': BCType.PRESSURE_OUTLET,
            'symmetry': BCType.SYMMETRY,
            'symmetryPlane': BCType.SYMMETRY,
            'cyclic': BCType.PERIODIC,
            'empty': BCType.EMPTY,
        }
        return mapping.get(of_type.lower(), BCType.UNKNOWN)

    def _parse_fields(self, zero_dir: Path, result: ConversionResult) -> int:
        """Parse field files from the 0/ directory.

        Reads ``internalField`` only (uniform and nonuniform List<...>). Files
        whose internal field cannot be read are recorded in the gap report
        rather than dropped. Returns the number of fields parsed.
        """
        parsed = 0
        for field_file in sorted(zero_dir.iterdir()):
            if not field_file.is_file() or field_file.name.startswith('.'):
                continue
            try:
                content = field_file.read_text()
            except OSError as exc:
                result.gap_report.unavailable(
                    "fields", field_file.name, f"Cannot read field file: {exc}"
                )
                continue

            match = re.search(r"internalField\s+([^;]+);", content, re.DOTALL)
            if not match:
                result.gap_report.unavailable(
                    "fields",
                    field_file.name,
                    "No internalField entry found (boundary-only or unsupported syntax)",
                )
                continue

            name = field_file.name
            parsed_field = self._parse_field_value(match.group(1))
            if parsed_field is None:
                result.gap_report.unsupported_nonblocking(
                    "fields",
                    name,
                    f"Unsupported internalField syntax: {' '.join(match.group(1).split())[:80]}",
                    "Re-export the field with a uniform or nonuniform List<scalar/vector> internalField",
                )
                continue

            if parsed_field.ndim == 1:
                result.scalar_fields.append((name, parsed_field))
            else:
                result.vec_fields.append((name, parsed_field))
            parsed += 1

        return parsed

    def _parse_field_value(self, text: str) -> np.ndarray | None:
        """Parse an OpenFOAM internalField value into a numpy array."""
        text = re.sub(r"/\*.*?\*/", " ", text, flags=re.DOTALL)
        text = re.sub(r"//.*$", "", text, flags=re.MULTILINE).strip()
        if not text:
            return None

        uniform = re.match(r"^uniform\s+(.+?)\s*$", text, re.DOTALL)
        if uniform:
            tokens = uniform.group(1).replace("(", " ").replace(")", " ").split()
            try:
                values = [float(t) for t in tokens]
            except ValueError:
                return None
            if not values:
                return None
            arr = np.array(values, dtype=np.float64)
            return arr

        nonuniform = re.match(r"^nonuniform\s+(.*)$", text, re.DOTALL)
        body = nonuniform.group(1) if nonuniform else text

        block = re.match(
            r"^(?:List<(\w+)>|(\w+))\s*(\d+)?\s*\((.*)\)\s*$", body, re.DOTALL
        )
        if not block:
            return None
        kind = block.group(1) or block.group(2)
        dims = _FIELD_DIMS.get(kind)
        if dims is None:
            return None

        entries = [e for e in re.split(r"[()]", block.group(4)) if e.strip()]
        values: list[float] = []
        for entry in entries:
            tokens = entry.replace(",", " ").split()
            try:
                values.extend(float(t) for t in tokens)
            except ValueError:
                return None
        if not values or len(values) % dims != 0:
            return None
        arr = np.array(values, dtype=np.float64)
        return arr if dims == 1 else arr.reshape(-1, dims)