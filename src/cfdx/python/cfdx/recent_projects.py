"""Bounded, deterministic recent-case preferences for Workbench.

This module never opens projects, deserializes executable content, or scans
folders. Only successfully loaded/saved project paths are recorded by Qt.
"""
from __future__ import annotations

from pathlib import Path

MAX_RECENTS = 8


def normalize_recent_projects(raw: object, *, require_existing: bool = True) -> tuple[str, ...]:
    """Deduplicate, validate and bound settings entries without globbing."""
    if not isinstance(raw, (list, tuple)):
        return ()
    selected: list[str] = []
    seen: set[str] = set()
    for value in raw[:64]:
        if not isinstance(value, str) or not 0 < len(value) <= 4096 or "\x00" in value:
            continue
        candidate = Path(value).expanduser()
        if not candidate.is_absolute() or candidate.suffix.lower() not in {".h5", ".cfdx"}:
            continue
        try:
            # Never offer symlinks as reusable entry points to a case.
            if candidate.is_symlink() or (require_existing and not candidate.is_file()):
                continue
            normalized = str(candidate.resolve(strict=require_existing))
        except (OSError, RuntimeError, ValueError):
            continue
        if normalized not in seen:
            seen.add(normalized)
            selected.append(normalized)
        if len(selected) >= MAX_RECENTS:
            break
    return tuple(selected)


def remember_recent_project(path: str | Path, previous: object) -> tuple[str, ...]:
    """Move a *confirmed* project to the front; ignore invalid targets."""
    candidate = normalize_recent_projects([str(Path(path))])
    if not candidate:
        return normalize_recent_projects(previous)
    return normalize_recent_projects([*candidate, *normalize_recent_projects(previous)])
