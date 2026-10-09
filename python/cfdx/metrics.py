"""Parsing of solver progress metrics from CFDX runner output."""
from __future__ import annotations

from dataclasses import dataclass
import re


@dataclass(frozen=True)
class SolverMetrics:
    """Latest metrics extracted from one solver output line."""

    iteration: int | None = None
    time: float | None = None
    cfl: float | None = None
    residuals: tuple[tuple[str, float], ...] = ()

    def residual(self, name: str) -> float | None:
        """Return a named residual, or None when it was not reported."""
        for key, value in self.residuals:
            if key == name:
                return value
        return None


class SolverMetricsParser:
    """Extract iteration, physical time, CFL and residuals from text lines.

    The parser is deliberately format-tolerant: labels may use '=' or ':'
    and numeric values may use either decimal or scientific notation. Unknown
    lines are ignored.
    """

    _NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    # Production solver completion reports "Iterations N" (plural), while
    # progress lines use "Iteration N". Accept both forms so the final
    # authoritative solver iteration is not lost by the monitor parser.
    _ITERATION = re.compile(
        rf"\b(?:iterations?|iter)\s*(?:=|:)?\s*(\d+)\b", re.I
    )
    _TIME = re.compile(rf"\b(?:time|physical\s+time)\s*(?:=|:)\s*({_NUMBER})\b", re.I)
    _CFL = re.compile(rf"\bCFL\s*(?:=|:)\s*({_NUMBER})\b", re.I)
    _RESIDUAL = re.compile(
        rf"\b(?:residual|res)\s*[(:]?\s*([A-Za-z0-9_.-]+)\s*[)]?\s*(?:=|:)\s*({_NUMBER})\b",
        re.I,
    )

    def parse(self, line: str) -> SolverMetrics | None:
        """Parse one line and return metrics when at least one is present."""
        iteration = self._match_int(self._ITERATION, line)
        time = self._match_float(self._TIME, line)
        cfl = self._match_float(self._CFL, line)

        residuals: list[tuple[str, float]] = []
        for match in self._RESIDUAL.finditer(line):
            residuals.append((match.group(1), float(match.group(2))))

        if iteration is None and time is None and cfl is None and not residuals:
            return None
        return SolverMetrics(iteration, time, cfl, tuple(residuals))

    @staticmethod
    def _match_int(pattern: re.Pattern[str], line: str) -> int | None:
        match = pattern.search(line)
        return int(match.group(1)) if match else None

    @staticmethod
    def _match_float(pattern: re.Pattern[str], line: str) -> float | None:
        match = pattern.search(line)
        return float(match.group(1)) if match else None
