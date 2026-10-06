"""Project/case lifecycle boundary shared by GUI, TUI and CLI.

A Project owns the canonical CFDX case artifact and its paired numerical DAT
restart artifact. It deliberately does not introduce a project manifest:
.cfdx.h5 remains setup/mesh state and .dat.h5 remains numerical state.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .case_io import (
    _paired_dat_path,
    read_case,
    read_case_with_dat,
    save_case,
    save_case_with_dat,
    validate_case_bundle,
)
from .session import CFDXSession


@dataclass(frozen=True)
class Project:
    """Canonical on-disk identity of one CFDX case.

    ``path`` is always the authoritative ``*.cfdx.h5`` setup artifact.
    The numerical restart, when present, is its sibling ``*.dat.h5``.
    """

    path: Path

    def __post_init__(self) -> None:
        path = Path(self.path)
        if not path.name.lower().endswith(".cfdx.h5"):
            raise ValueError("CFDX project case must use the canonical .cfdx.h5 extension")
        object.__setattr__(self, "path", path)

    @property
    def dat_path(self) -> Path:
        """Return the canonical sibling numerical checkpoint path."""
        return _paired_dat_path(self.path)

    @property
    def name(self) -> str:
        """Return the project stem without CFDX artifact suffixes."""
        return self.path.name[: -len(".cfdx.h5")]

    def save(self, session: CFDXSession) -> Path:
        """Persist setup/mesh state only; never consumes or writes DAT state."""
        return save_case(session, self.path)

    def save_with_dat(self, session: CFDXSession, source_dat: Path) -> tuple[Path, Path]:
        """Persist the case and canonicalize an explicit numerical checkpoint."""
        return save_case_with_dat(session, self.path, Path(source_dat))

    def load(self) -> CFDXSession:
        """Load setup/mesh state without loading numerical restart state."""
        return read_case(self.path)

    def load_with_dat(self, dat_path: Path | None = None) -> tuple[CFDXSession, Path]:
        """Load setup and explicitly restore the paired numerical state."""
        return read_case_with_dat(self.path, dat_path)

    def validate(self) -> dict[str, bool]:
        """Validate the self-contained CFDX case artifact."""
        return validate_case_bundle(self.path)

    def exists(self) -> bool:
        return self.path.is_file()

    def has_dat(self) -> bool:
        return self.dat_path.is_file()
