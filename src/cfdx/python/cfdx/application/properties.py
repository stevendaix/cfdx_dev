"""Contextual property contracts for the Workbench Properties panel.

The same resolver serves every frontend: workflow selections resolve to
editable typed properties owned by the Case, while mesh selections
(regions, patches, faces) resolve to read-only typed facts owned by the
C++/HDF5 mesh topology.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..case import Case
from ..mesh_model import MeshCatalog, MeshPatch, MeshRegion
from ..selection import FACE_KIND, REGION_KIND, SelectionTarget, parse_selection_id
from ..setup_model import ChangeImpact, Parameter, ParameterType


@dataclass(frozen=True)
class PropertyState:
    key: str
    label: str
    value: Any
    kind: ParameterType
    unit: str | None = None
    editable: bool = True
    impact: ChangeImpact = ChangeImpact.HOT
    message: str | None = None


def properties_for_selection(
    case: Case, stable_id: str | None, catalog: MeshCatalog | None = None
) -> tuple[PropertyState, ...]:
    """Resolve a stable navigation ID to typed application properties.

    Mesh selection ids (``region:``/``patch:``) are delegated to the mesh
    resolver and raise when they reference objects missing from the
    catalog; workflow ids never raise.
    """
    if stable_id is None:
        return ()
    if stable_id.startswith("region:") or stable_id.startswith("patch:"):
        if catalog is None:
            return ()
        return mesh_selection_properties(catalog, stable_id)
    if stable_id == "numerics":
        value = case.numerics.get("cfl", 1.0)
        return (
            PropertyState("numerics.cfl", "CFL", value, ParameterType.REAL),
        )
    if stable_id == "solver":
        return (
            PropertyState(
                "execution.solver",
                "Solver executable",
                case.execution.solver or "",
                ParameterType.STRING,
                impact=ChangeImpact.REQUIRES_RESTART,
            ),
        )
    if stable_id == "physics":
        return tuple(
            PropertyState(
                f"physics.{name}.enabled",
                name.title(),
                bool(value.get("enabled", True)),
                ParameterType.BOOL,
                impact=ChangeImpact.REQUIRES_RESTART,
            )
            for name, value in case.physics.items()
            if isinstance(value, dict)
        )
    if stable_id == "materials":
        return tuple(
            PropertyState(
                f"materials.{name}",
                name.title(),
                value,
                ParameterType.STRING,
                editable=False,
            )
            for name, value in case.materials.items()
        )
    if stable_id == "boundaries":
        return tuple(
            PropertyState(
                f"boundaries.{name}.type",
                name,
                value.get("type", ""),
                ParameterType.STRING,
                editable=False,
            )
            for name, value in case.boundaries.items()
            if isinstance(value, dict)
        )
    return ()


def mesh_selection_properties(
    catalog: MeshCatalog, stable_id: str
) -> tuple[PropertyState, ...]:
    """Resolve one canonical mesh selection id to read-only typed rows.

    Mesh topology is owned by the C++/HDF5 layer: the inspector may
    describe it, never offer to edit it, so every row is read-only.
    """
    target = parse_selection_id(stable_id)
    if target.kind == REGION_KIND:
        return _region_properties(_region_in_catalog(catalog, target.name))
    patch = _patch_in_catalog(catalog, target)
    if target.kind == FACE_KIND:
        if target.face_index is None or target.face_index >= patch.face_count:
            raise ValueError(f"mesh face out of range: {target.stable_id}")
        return (
            _read_only("mesh.object", "Object type", "Face", ParameterType.STRING),
            _read_only(
                "mesh.face", "Face index", target.face_index, ParameterType.INTEGER
            ),
            _read_only("mesh.patch", "Patch", patch.name, ParameterType.STRING),
            _read_only("mesh.region", "Region", patch.region, ParameterType.STRING),
        )
    return _patch_properties(patch)


def set_property(case: Case, key: str, value: Any) -> ChangeImpact:
    """Apply an explicitly routed property to the authoritative Case."""
    parts = key.split(".")
    if key == "numerics.cfl":
        parameter = Parameter("cfl", value, ParameterType.REAL)
        parameter.validate()
        case.numerics["cfl"] = value
        return ChangeImpact.HOT
    if key == "execution.solver":
        if not isinstance(value, str):
            raise ValueError("solver executable must be a string")
        case.execution.solver = value or None
        return ChangeImpact.REQUIRES_RESTART
    if len(parts) == 3 and parts[0] == "physics" and parts[2] == "enabled":
        if not isinstance(value, bool):
            raise ValueError("physics enabled must be boolean")
        model = case.physics.setdefault(parts[1], {})
        if not isinstance(model, dict):
            raise ValueError(f"physics.{parts[1]} must be a mapping")
        model["enabled"] = value
        return ChangeImpact.REQUIRES_RESTART
    raise KeyError(f"property is not editable: {key}")


def _read_only(
    key: str, label: str, value: Any, kind: ParameterType
) -> PropertyState:
    """Build one read-only inspector row for a mesh fact."""
    return PropertyState(key, label, value, kind, editable=False)


def _region_properties(region: MeshRegion) -> tuple[PropertyState, ...]:
    return (
        _read_only("mesh.object", "Object type", "Region", ParameterType.STRING),
        _read_only("mesh.cells", "Cells", region.n_cells, ParameterType.INTEGER),
        _read_only("mesh.faces", "Faces", region.n_faces, ParameterType.INTEGER),
        _read_only("mesh.points", "Points", region.n_points, ParameterType.INTEGER),
        _read_only(
            "mesh.patches", "Patches", len(region.patches), ParameterType.INTEGER
        ),
    )


def _patch_properties(patch: MeshPatch) -> tuple[PropertyState, ...]:
    return (
        _read_only("mesh.object", "Object type", "Patch", ParameterType.STRING),
        _read_only("mesh.patch_type", "Patch type", patch.type, ParameterType.STRING),
        _read_only("mesh.faces", "Faces", patch.face_count, ParameterType.INTEGER),
        _read_only("mesh.region", "Region", patch.region, ParameterType.STRING),
    )


def _region_in_catalog(catalog: MeshCatalog, name: str) -> MeshRegion:
    for region in catalog.regions:
        if region.name == name:
            return region
    raise KeyError(f"mesh region not found: {name}")


def _patch_in_catalog(
    catalog: MeshCatalog, target: SelectionTarget
) -> MeshPatch:
    region = target.region or "default"
    for patch in catalog.patches:
        if patch.name == target.name and patch.region == region:
            return patch
    raise KeyError(f"mesh patch not found: {target.stable_id}")
