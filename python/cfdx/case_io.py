"""CFDX case/restart artifact persistence.

The HDF5 case is the persistent source of truth. A solver DAT file is kept
as a paired numerical-restart artifact when explicitly saved.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import h5py

from .case import Case, ExecutionConfig
from .dat_io import read_dat_restart, write_dat_hdf5
from .session import CFDXSession


_FORMAT = "CFDX"
_SCHEMA_VERSION = 1
_CASE_DATASET = "case/config"
_CHECKPOINT_GROUP = "runtime/checkpoint"
_RESTART_GROUP = "runtime/restart"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_case_hash(case_data: dict, mesh_hash: str | None = None) -> str:
    """Hash normalized case configuration plus the persisted mesh identity."""
    payload = {"case": case_data, "mesh_hash": mesh_hash}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_path(path: Path) -> Path:
    path = Path(path)
    if path.suffix.lower() != ".h5":
        raise ValueError("CFDX case must use an .h5 file")
    return path


def _paired_dat_path(case_path: Path) -> Path:
    """Return the canonical DAT sibling for case.cfdx.h5."""
    name = case_path.name
    if name.lower().endswith(".cfdx.h5"):
        name = name[:-len(".cfdx.h5")]
    else:
        name = case_path.stem
    return case_path.with_name(f"{name}.dat")


def save_case(session: CFDXSession, path: Path) -> Path:
    """Save the case configuration and application checkpoint to HDF5."""
    path = _validate_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    case_data = session.case.as_dict()

    # Append to an existing CFDX HDF5 artifact so a mesh/topology written by
    # the native C++ HDF5 layer is never destroyed by a case save.
    with h5py.File(path, "a") as h5:
        h5.attrs["format"] = _FORMAT
        h5.attrs["schema_version"] = _SCHEMA_VERSION
        h5.attrs["case_revision"] = session.case_revision
        h5.attrs["mesh_revision"] = session.mesh_revision
        h5.attrs["physics_revision"] = session.physics_revision
        h5.attrs["numerics_revision"] = session.numerics_revision
        now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        if "creation_date" not in h5.attrs:
            h5.attrs["creation_date"] = now
        h5.attrs["modification_date"] = now
        h5.attrs["dimension"] = int(h5.attrs.get("dimension", 3))
        h5.attrs["precision"] = str(h5.attrs.get("precision", "float64"))
        h5.attrs["endian"] = str(h5.attrs.get("endian", "little"))

        case_group = h5.require_group("case")
        if "config" in case_group:
            del case_group["config"]
        case_group.create_dataset(
            "config",
            data=json.dumps(case_data, sort_keys=True, separators=(",", ":")),
        )

        mesh_hash = h5.attrs.get("mesh_hash")
        if isinstance(mesh_hash, bytes):
            mesh_hash = mesh_hash.decode("utf-8")
        h5.attrs["case_hash"] = _canonical_case_hash(case_data, mesh_hash)
        if "geometry_hash" in h5.attrs:
            h5.attrs["geometry_hash"] = str(h5.attrs["geometry_hash"])

        checkpoint = h5.require_group(_CHECKPOINT_GROUP)
        checkpoint.attrs["iteration"] = session.iteration
        checkpoint.attrs["time"] = session.time

    return path


def save_case_with_dat(
    session: CFDXSession, path: Path, dat_path: Path
) -> tuple[Path, Path]:
    """Save HDF5 case state and copy the solver DAT beside it."""
    dat_path = Path(dat_path)
    if not dat_path.is_file():
        raise FileNotFoundError(dat_path)

    case_path = save_case(session, path)
    target_dat = _paired_dat_path(case_path)
    restart = read_dat_restart(dat_path)
    write_dat_hdf5(target_dat, restart)

    with h5py.File(case_path, "r+") as h5:
        restart = h5.require_group(_RESTART_GROUP)
        restart.attrs["dat_filename"] = target_dat.name
        restart.attrs["dat_sha256"] = _sha256(target_dat)
        restart.attrs["dat_size"] = target_dat.stat().st_size

    return case_path, target_dat


def _read_session(path: Path) -> CFDXSession:
    path = _validate_path(path)
    if not path.is_file():
        raise FileNotFoundError(path)

    with h5py.File(path, "r") as h5:
        file_format = h5.attrs.get("format")
        if isinstance(file_format, bytes):
            file_format = file_format.decode("utf-8")
        if file_format != _FORMAT:
            raise ValueError("not a CFDX HDF5 case")
        if int(h5.attrs.get("schema_version", -1)) != _SCHEMA_VERSION:
            raise ValueError("unsupported CFDX HDF5 case schema")

        try:
            raw = h5[_CASE_DATASET][()]
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            data = json.loads(raw)
            execution = ExecutionConfig(**data["execution"])
            case = Case(
                name=data["name"],
                physics=dict(data["physics"]),
                numerics=dict(data["numerics"]),
                boundaries=dict(data["boundaries"]),
                materials=dict(data.get("materials", {})),
                execution=execution,
            )
            checkpoint = h5[_CHECKPOINT_GROUP].attrs
            session = CFDXSession(
                case=case,
                case_revision=int(h5.attrs["case_revision"]),
                mesh_revision=int(h5.attrs["mesh_revision"]),
                physics_revision=int(h5.attrs["physics_revision"]),
                numerics_revision=int(h5.attrs["numerics_revision"]),
                iteration=int(checkpoint["iteration"]),
                time=float(checkpoint["time"]),
            )
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("invalid CFDX HDF5 case content") from exc

    return session


def validate_case_bundle(path: Path) -> dict[str, bool]:
    """Validate one self-contained CFDX HDF5 case artifact."""
    case_path = _validate_path(path)
    if not case_path.is_file():
        raise FileNotFoundError(case_path)
    required_mesh = ("points", "face_vertices", "face_offsets", "owner", "neighbour", "cell_faces", "cell_offsets")
    with h5py.File(case_path, "r") as h5:
        if h5.attrs.get("format") != _FORMAT:
            raise ValueError("not a CFDX HDF5 case")
        if int(h5.attrs.get("schema_version", -1)) != _SCHEMA_VERSION:
            raise ValueError("unsupported CFDX HDF5 case schema")
        if _CASE_DATASET not in h5:
            raise ValueError("case configuration is missing")
        if not all(name in h5 for name in required_mesh):
            raise ValueError("mesh data is missing from the self-contained case")
        if "case_hash" not in h5.attrs:
            raise ValueError("case integrity hash is missing")
        mesh_hash = h5.attrs.get("mesh_hash")
        if isinstance(mesh_hash, bytes):
            mesh_hash = mesh_hash.decode("utf-8")
        expected_case_hash = _canonical_case_hash(data, mesh_hash)
        recorded_case_hash = h5.attrs["case_hash"]
        if isinstance(recorded_case_hash, bytes):
            recorded_case_hash = recorded_case_hash.decode("utf-8")
        if str(recorded_case_hash) != expected_case_hash:
            raise ValueError("case integrity hash mismatch")
        if "fields" not in h5 or "values" not in h5["fields"]:
            raise ValueError("fields data is missing from the self-contained case")
        raw = h5[_CASE_DATASET][()]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        data = json.loads(raw)
        if not isinstance(data.get("physics"), dict):
            raise ValueError("physics configuration is missing")
        if not isinstance(data.get("numerics"), dict):
            raise ValueError("numerics configuration is missing")
        runtime_separated = "runtime" in h5 and "case" in h5 and "config" not in h5["runtime"]
    return {"mesh": True, "fields": True, "physics": True, "numerics": True, "runtime_separated": runtime_separated}


def read_case(path: Path) -> CFDXSession:
    """Load only the HDF5 case; no solver restart artifact is consumed."""
    return _read_session(path)


def read_case_with_dat(
    path: Path, dat_path: Path | None = None
) -> tuple[CFDXSession, Path]:
    """Load a case and validate its paired solver DAT artifact."""
    case_path = _validate_path(path)
    session = _read_session(case_path)

    with h5py.File(case_path, "r") as h5:
        if _RESTART_GROUP not in h5:
            raise ValueError("case does not contain a paired DAT restart")
        restart = h5[_RESTART_GROUP].attrs
        recorded_name = str(restart["dat_filename"])
        recorded_hash = str(restart["dat_sha256"])

    candidate = Path(dat_path) if dat_path is not None else case_path.with_name(recorded_name)
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    if _sha256(candidate) != recorded_hash:
        raise ValueError("DAT restart is incompatible with the CFDX case")

    return session, candidate
