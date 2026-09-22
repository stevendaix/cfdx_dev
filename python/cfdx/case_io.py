"""CFDX case/restart artifact persistence.

The HDF5 case is the persistent source of truth. A solver DAT file is kept
as a paired numerical-restart artifact when explicitly saved.
"""

from __future__ import annotations

import hashlib
import json
import shutil
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

        case_group = h5.require_group("case")
        if "config" in case_group:
            del case_group["config"]
        case_group.create_dataset(
            "config",
            data=json.dumps(case_data, sort_keys=True, separators=(",", ":")),
        )

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
