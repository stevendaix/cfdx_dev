"""Stable selection identity scheme shared by the GUI browser and 3D view.

The mesh browser and the PyVistaQt view navigate the same mesh objects, so
both consume one canonical id scheme: ``region:<name>`` for regions,
``[region:<r>/]patch:<name>`` for patches and ``<patch id>/face:<k>`` for
zero-based faces.  The default region carries no prefix, ids are strict,
and anything the parser rejects can never reach the renderer.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

PATCH_KIND = "patch"
FACE_KIND = "face"
REGION_KIND = "region"

_INVALID_ID = "invalid selection stable id"

_FACE_PATTERN = (
    r"^(?:(?P<region_kind>region):(?P<region>[^/:]+)/)?"
    r"(?P<patch_kind>patch):(?P<patch>[^/:]+)/"
    r"(?P<face_kind>face):(?P<face_index>(0|[1-9][0-9]*))$"
)
_PATCH_PATTERN = (
    r"^(?:(?P<region_kind>region):(?P<region>[^/:]+)/)?"
    r"(?P<patch_kind>patch):(?P<patch>[^/:]+)$"
)
_REGION_PATTERN = r"^region:(?P<region>[^/:]+)$"

_FACE_RE = re.compile(_FACE_PATTERN)
_PATCH_RE = re.compile(_PATCH_PATTERN)
_REGION_RE = re.compile(_REGION_PATTERN)


@dataclass(frozen=True)
class SelectionTarget:
    """One canonical mesh selection: a region, a patch or a single face."""

    kind: str
    stable_id: str
    name: str
    region: str | None = None
    face_index: int | None = None


def _canonical_name(value: str) -> str:
    """Reject names that would break the canonical id syntax."""
    if not value or value != value.strip() or "/" in value or ":" in value:
        raise ValueError(_INVALID_ID)
    return value


def parse_selection_id(stable_id: str) -> SelectionTarget:
    """Parse one canonical selection id or raise ``ValueError``."""
    match = _FACE_RE.match(stable_id)
    if match is not None:
        region = match.group("region")
        if region is not None:
            _canonical_name(region)
        return SelectionTarget(
            FACE_KIND,
            stable_id,
            _canonical_name(match.group("patch")),
            region=region,
            face_index=int(match.group("face_index")),
        )
    match = _PATCH_RE.match(stable_id)
    if match is not None:
        region = match.group("region")
        if region is not None:
            _canonical_name(region)
        return SelectionTarget(
            PATCH_KIND,
            stable_id,
            _canonical_name(match.group("patch")),
            region=region,
        )
    match = _REGION_RE.match(stable_id)
    if match is not None:
        return SelectionTarget(
            REGION_KIND, stable_id, _canonical_name(match.group("region"))
        )
    raise ValueError(_INVALID_ID)


def patch_selection_id(name: str, region: str | None = None) -> str:
    """Build the canonical patch id; the default region adds no prefix."""
    canonical = _canonical_name(name)
    if region is None or region == "default":
        return f"patch:{canonical}"
    return f"region:{_canonical_name(region)}/patch:{canonical}"


def region_selection_id(name: str) -> str:
    """Build the canonical region id."""
    return f"region:{_canonical_name(name)}"


def face_selection_id(patch_id: str, face_index: int) -> str:
    """Build the canonical zero-based face id for a loaded patch id."""
    if parse_selection_id(patch_id).kind != PATCH_KIND:
        raise ValueError(_INVALID_ID)
    if isinstance(face_index, bool) or not isinstance(face_index, int):
        raise ValueError(_INVALID_ID)
    if face_index < 0:
        raise ValueError(_INVALID_ID)
    return f"{patch_id}/face:{face_index}"


def selection_label(target: SelectionTarget) -> str:
    """Return a short human-readable label for status bars."""
    if target.kind == REGION_KIND:
        return f"Region {target.name}"
    if target.kind == FACE_KIND:
        patch_part = f"Face {target.face_index} on patch {target.name}"
        if target.region is None:
            return patch_part
        return f"{patch_part} — region {target.region}"
    if target.kind == PATCH_KIND:
        if target.region is None:
            return f"Patch {target.name}"
        return f"Patch {target.name} — region {target.region}"
    raise ValueError(_INVALID_ID)
