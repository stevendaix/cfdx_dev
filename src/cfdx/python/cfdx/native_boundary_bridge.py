"""Explicit, lossless Workbench field-boundary to native CFDX C++ CaseSetup bridge.

Only a native case with existing case_setup_json can be the source. Preserve
all numerical, turbulence, and material settings; do not invent them from a
GUI form. A new output case is created transactionally (never in-place).
"""

from __future__ import annotations

import json
import math
import os
import shutil
import tempfile
from collections.abc import Mapping
from pathlib import Path

import h5py

from .boundary_editor import condition_for_field

_ALLOWED = {"velocity", "pressure", "temperature", "k", "omega"}
_KIND = {"inlet", "outlet", "wall", "symmetry", "pressure_outlet"}


class NativeBoundaryError(ValueError):
    """Editor data cannot be applied to the native C++ solver safely."""


def _finite(value, name):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise NativeBoundaryError(f"{name} must be finite")
    return float(value)


def translate_patch(
    name: str, editor: Mapping[str, object], original: Mapping[str, object]
) -> dict:
    """Translate all declared fields; preserve unrelated native BC metadata."""
    kind = editor.get("type")
    if kind not in _KIND:
        raise NativeBoundaryError(f"unsupported native patch kind for {name}: {kind}")
    if kind != original.get("type") and not (
        kind == "outlet" and original.get("type") == "pressure_outlet"
    ):
        raise NativeBoundaryError(
            f"{name}: patch kind differs from native case; do not reinterpret geometry"
        )
    fields = editor.get("fields", [])
    if (
        not isinstance(fields, list)
        or len(fields) != len(set(fields))
        or any(field not in _ALLOWED for field in fields)
    ):
        raise NativeBoundaryError(
            f"{name}: duplicate, unsupported or malformed field list"
        )
    if kind == "symmetry" and fields:
        raise NativeBoundaryError(
            f"{name}: symmetry cannot declare explicit fields in this bridge"
        )
    result = dict(original)
    # Preserve unmodified SST conditions and case provenance. Do not project a
    # velocity condition onto pressure or temperature via the legacy value_type.
    native_fields = dict(result.get("field_conditions", {}))
    sst_conditions = {
        item["field"]: dict(item) for item in result.get("sst_scalar_conditions", [])
    }
    for field in fields:
        info = condition_for_field(editor, field)
        mode = info.get("velocity_type" if field == "velocity" else "scalar_type")
        if field == "velocity":
            if mode not in {"FIXED_VALUE", "ZERO_GRADIENT"}:
                raise NativeBoundaryError(f"{name}: velocity mode is not supported")
            velocity = info.get("velocity_value")
            if not isinstance(velocity, (tuple, list)) or len(velocity) != 3:
                raise NativeBoundaryError(f"{name}: velocity requires three components")
            values = [_finite(v, f"{name}.velocity") for v in velocity]
            native_fields[field] = {"type": mode, "value": values, "gradient": 0.0}
        elif field in ("pressure", "temperature"):
            if mode not in {"FIXED_VALUE", "ZERO_GRADIENT", "FIXED_GRADIENT"}:
                raise NativeBoundaryError(f"{name}: {field} mode unsupported")
            value = _finite(info.get("value", 0), f"{name}.{field}.value")
            gradient = _finite(info.get("gradient", 0), f"{name}.{field}.gradient")
            if field == "temperature" and mode == "FIXED_VALUE" and value <= 0:
                raise NativeBoundaryError(
                    "absolute temperature must be positive kelvin"
                )
            native_fields[field] = {
                "type": mode,
                "value": [value],
                "gradient": gradient,
            }
        elif field in ("k", "omega"):
            value = _finite(info.get("value", 0), f"{name}.{field}.value")
            if kind == "wall":
                if (
                    field == "k"
                    and mode in ("WALL_K_ZERO", "FIXED_VALUE")
                    and value == 0
                ):
                    native_type = "wall_k_zero"
                elif field == "omega" and mode == "WALL_OMEGA_SST":
                    native_type = "wall_omega_sst"
                else:
                    raise NativeBoundaryError(
                        f"{name}: wall {field} requires a certified wall condition"
                    )
            elif mode in {"FIXED_VALUE", "ZERO_GRADIENT"}:
                if mode == "FIXED_VALUE" and value <= 0:
                    raise NativeBoundaryError(
                        f"{field} fixed boundary must be positive"
                    )
                native_type = mode.lower()
            else:
                raise NativeBoundaryError(f"{name}: unsupported SST mode for {field}")
            sst_conditions[field] = {
                "field": field,
                "type": native_type,
                "value": value,
            }
    result["field_conditions"] = native_fields
    result["sst_scalar_conditions"] = [
        sst_conditions[key] for key in sorted(sst_conditions)
    ]
    return result


def build_native_boundary_payload(
    existing: dict, editor_boundaries: Mapping[str, object]
) -> dict:
    """Update only known patch boundaries in the existing native C++ JSON."""
    if not isinstance(existing, dict) or not isinstance(
        existing.get("boundary_conditions"), list
    ):
        raise NativeBoundaryError("native CaseSetup boundary_conditions missing")
    original = existing["boundary_conditions"]
    names = [entry.get("patch_name") for entry in original if isinstance(entry, dict)]
    if len(names) != len(set(names)) or any(
        not isinstance(name, str) or not name for name in names
    ):
        raise NativeBoundaryError("invalid or duplicated native patch names")
    if set(editor_boundaries) != set(names):
        raise NativeBoundaryError(
            "editor and native mesh patches differ; all native patches must be configured"
        )
    result = dict(existing)
    result["boundary_conditions"] = [
        translate_patch(name, editor_boundaries[name], entry)
        for name, entry in zip(names, original)
    ]
    return result


def write_native_boundary_case(
    source: str | Path, destination: str | Path, editor_boundaries: Mapping[str, object]
) -> Path:
    """Create a separate, solver-readable HDF5 artifact; leave source untouched.

    This writes the native C++ root attribute, not /case/config. A caller must
    clearly label the result an exported *solver case*, not an editor save.
    """
    source, destination = Path(source).resolve(), Path(destination).absolute()
    if not source.is_file() or source == destination or destination.is_symlink():
        raise NativeBoundaryError(
            "source must exist and output must be a different, non-symlink path"
        )
    if not destination.name.endswith(".cfdx.h5"):
        raise NativeBoundaryError("native solver case output must end in .cfdx.h5")
    if destination.exists():
        raise NativeBoundaryError(
            "destination exists; refuse to overwrite another case"
        )
    with h5py.File(source, "r") as h5:
        raw = h5.attrs.get("case_setup_json")
        if not isinstance(raw, (str, bytes)):
            raise NativeBoundaryError(
                "source has no native C++ case_setup_json; cannot synthesize physics"
            )
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        original = json.loads(raw)
    payload = build_native_boundary_payload(original, editor_boundaries)
    data = json.dumps(payload, separators=(",", ":"), sort_keys=True, allow_nan=False)
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle, temp = tempfile.mkstemp(
        prefix=".cfdx-native-bcs-", suffix=".cfdx.h5", dir=destination.parent
    )
    os.close(handle)
    try:
        shutil.copyfile(source, temp)  # stream copy, not load whole mesh into RAM
        with h5py.File(temp, "r+") as h5:
            h5.attrs["case_setup_json"] = data
            h5.attrs["workbench_native_boundary_bridge_version"] = "1"
        os.replace(temp, destination)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)
    return destination
