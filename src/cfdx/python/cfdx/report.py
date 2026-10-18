"""Deterministic, exportable text reports for CFDX runs."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Report:
    title: str
    sections: Mapping[str, str]

    def to_markdown(self) -> str:
        lines = [f"# {self.title}", ""]
        for name, body in self.sections.items():
            lines.extend([f"## {name}", "", body.strip(), ""])
        return "\n".join(lines)

    def write_markdown(self, path: Path) -> None:
        path.write_text(self.to_markdown(), encoding="utf-8")
