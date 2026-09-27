"""OpenFOAM case directory adapter for CFDX I/O.

Parses OpenFOAM case directory structure:
  - constant/polyMesh/points, faces, owner, neighbour, boundary
  - 0/ field files (U, p, T, etc.)

Corresponds to the C++ header:
  src/cfdx/io/openfoam/openfoam_importer.h
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union

import numpy as np

from cfdx.io.interfaces import (
    SolverAdapter,
    ConversionResult,
)
from cfdx.io.schema import (
    BCType,
    BCValueType,
    CaseSetup,
    MaterialSpec,
    BoundarySpec,
    SourceInfo,
)
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.gap_analysis import GapAnalysis


_FIELD_DIMS = {"scalar": 1, "vector": 3, "symmTensor": 6, "tensor": 9}


def _extract_foam_list(content: str):
    """Extract the first parenthesised list block from a FoamFile ascii file.

    Returns ``(lines, count)`` where ``lines`` are the comment-stripped,
    non-empty lines inside the block and ``count`` is the integer written just
    before the opening parenthesis (``None`` when absent).  Returns
    ``(None, None)`` when no list block is present.
    """
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

    def __init__(self, rules: Optional[MappingRules] = None) -> None:
        self.rules = rules or MappingRules.load_default()
        self.setup = CaseSetup(source=SourceInfo(solver="OpenFOAM", format="openfoam"))

    def convert(self, case_path: Union[str, Path], result: ConversionResult) -> bool:
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
        points, cells = self._parse_mesh(poly_mesh)
        if points is None or cells is None:
            result.gap_report.unsupported_blocking(
                "mesh", "openfoam_polymesh_parse", "Failed to parse OpenFOAM mesh"
            )
            return False

        result.mesh = {"points": points, "cells": cells}
        result.gap_report.supported("mesh", "polyMesh topology")

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

        result.setup = self.setup
        result.source.case_path = str(path)
        result.source.solver = "OpenFOAM"
        result.source.format = "openfoam"

        return result.success

    def _parse_mesh(self, poly_mesh: Path) -> tuple[Optional[np.ndarray], Optional[dict]]:
        """Parse OpenFOAM polyMesh files."""
        try:
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

        except Exception as e:
            return None, None

    def _parse_points(self, points_file: Path) -> Optional[np.ndarray]:
        """Parse OpenFOAM points file."""
        try:
            content = points_file.read_text()

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
        except Exception:
            return None

            data_lines = match.group(1).strip().split('\n')
            # First line might be the count
            first_line = data_lines[0].strip()
            if first_line.isdigit():
                data_lines = data_lines[1:]

            points = []
            for line in data_lines:
                line = line.strip().rstrip(';')
                if line.startswith('(') and line.endswith(')'):
                    coords = line[1:-1].split()
                    if len(coords) >= 3:
                        points.append([float(coords[0]), float(coords[1]), float(coords[2])])

            return np.array(points, dtype=np.float64)
        except Exception:
            return None

    def _parse_faces(self, faces_file: Path) -> Optional[list]:
        """Parse OpenFOAM faces file.

        Handles both the compact ``4(0 1 2 3)`` form written by OpenFOAM and
        the parenthesised ``(0 1 2 3)`` form.
        """
        try:
            content = faces_file.read_text()

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
        except Exception:
            return None

    def _parse_label_list(self, label_file: Path) -> Optional[list]:
        """Parse OpenFOAM label list (owner, neighbour)."""
        try:
            content = label_file.read_text()

            data_lines, count = _extract_foam_list(content)
            if data_lines is None:
                return None

            labels = []
            for line in data_lines:
                token = line.strip().rstrip(';')
                if not token or token.startswith("/*") or token.startswith("//"):
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
        except Exception:
            return None

    def _parse_boundary_file(self, boundary_file: Path) -> list:
        """Parse OpenFOAM boundary file to get patch info."""
        patches = []
        try:
            with open(boundary_file, 'r') as f:
                content = f.read()

            # Parse each patch block
            patch_pattern = re.compile(
                r'(\w+)\s*\{[^}]*type\s+(\w+)[^}]*nFaces\s+(\d+)[^}]*startFace\s+(\d+)',
                re.DOTALL
            )
            for match in patch_pattern.finditer(content):
                name = match.group(1)
                patch_type = match.group(2)
                n_faces = int(match.group(3))
                start_face = int(match.group(4))
                patches.append({
                    'name': name,
                    'type': patch_type,
                    'n_faces': n_faces,
                    'start_face': start_face
                })
        except Exception:
            pass
        return patches

    def _build_cells_from_faces(self, faces, owner, neighbour, patches):
        """Build cell connectivity from face data."""
        # This is a simplified version - in reality we'd need to reconstruct
        # the full cell topology. For now, return a minimal structure.
        n_faces = len(faces)
        n_cells = max(owner) + 1 if owner else 0
        
        # Create a simple cell->faces mapping (each cell gets its owned faces)
        cell_faces = {i: [] for i in range(n_cells)}
        for face_idx, cell_idx in enumerate(owner):
            if cell_idx < n_cells:
                cell_faces[cell_idx].append(face_idx)

        # For meshio, we'll use a simpler representation - just return the faces
        # and let meshio handle the cell reconstruction
        # Return a simple triangle representation of faces for now
        cells = {}
        tri_cells = []
        for face in faces:
            if len(face) == 3:
                tri_cells.append(np.array(face, dtype=np.int64))
            elif len(face) == 4:
                # Split quad into two triangles
                tri_cells.append(np.array([face[0], face[1], face[2]], dtype=np.int64))
                tri_cells.append(np.array([face[0], face[2], face[3]], dtype=np.int64))

        if tri_cells:
            cells["triangle"] = np.array(tri_cells, dtype=np.int64)

        return cells

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

    def _parse_field_value(self, text: str) -> Optional[np.ndarray]:
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
            return arr if arr.size == 1 else arr

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