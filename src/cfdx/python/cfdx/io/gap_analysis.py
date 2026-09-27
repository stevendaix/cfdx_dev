"""Gap Analysis report engine — generates Markdown + JSON reports.

Corresponds to the C++ header:
  src/cfdx_io/gap_analysis.h/.cpp
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Optional

from cfdx.io.schema import Severity


@dataclass
class Finding:
    severity: Severity
    category: str
    feature: str
    detail: str = ""
    suggestion: str = ""


class GapAnalysis:
    """Accumulates findings from a conversion/import operation.

    Every unsupported or approximated feature must be recorded here.
    No silent substitutions/fallbacks.
    """

    def __init__(self) -> None:
        self._findings: list[Finding] = []

    def add(
        self,
        severity: Severity,
        category: str,
        feature: str,
        detail: str = "",
        suggestion: str = "",
    ) -> "GapAnalysis":
        self._findings.append(
            Finding(severity, category, feature, detail, suggestion)
        )
        return self

    def supported(self, category: str, feature: str, detail: str = "") -> "GapAnalysis":
        return self.add(Severity.SUPPORTED, category, feature, detail)

    def approximated(self, category: str, feature: str, detail: str, suggestion: str = "") -> "GapAnalysis":
        return self.add(Severity.APPROXIMATED, category, feature, detail, suggestion)

    def unsupported_nonblocking(self, category: str, feature: str, detail: str, suggestion: str = "") -> "GapAnalysis":
        return self.add(Severity.UNSUPPORTED_NONBLOCK, category, feature, detail, suggestion)

    def unsupported_blocking(self, category: str, feature: str, detail: str, suggestion: str = "") -> "GapAnalysis":
        return self.add(Severity.UNSUPPORTED_BLOCK, category, feature, detail, suggestion)

    def unavailable(self, category: str, feature: str, detail: str) -> "GapAnalysis":
        return self.add(Severity.UNAVAILABLE, category, feature, detail)

    @property
    def findings(self) -> list[Finding]:
        return list(self._findings)

    def has_blocking(self) -> bool:
        return any(f.severity == Severity.UNSUPPORTED_BLOCK for f in self._findings)

    def n_supported(self) -> int:
        return sum(1 for f in self._findings if f.severity in (Severity.SUPPORTED, Severity.APPROXIMATED))

    def n_approximated(self) -> int:
        return sum(1 for f in self._findings if f.severity == Severity.APPROXIMATED)

    def n_unsupported_nonblocking(self) -> int:
        return sum(1 for f in self._findings if f.severity == Severity.UNSUPPORTED_NONBLOCK)

    def n_unsupported_blocking(self) -> int:
        return sum(1 for f in self._findings if f.severity == Severity.UNSUPPORTED_BLOCK)

    def n_unavailable(self) -> int:
        return sum(1 for f in self._findings if f.severity == Severity.UNAVAILABLE)


class GapAnalysisReport:
    """Generates Markdown and JSON reports from a GapAnalysis."""

    def __init__(self, gap: GapAnalysis) -> None:
        self.gap = gap

    def to_json(self) -> str:
        data = {
            "has_blocking": self.gap.has_blocking(),
            "summary": {
                "supported": self.gap.n_supported(),
                "approximated": self.gap.n_approximated(),
                "unsupported_nonblocking": self.gap.n_unsupported_nonblocking(),
                "unsupported_blocking": self.gap.n_unsupported_blocking(),
                "unavailable": self.gap.n_unavailable(),
                "total": len(self.gap.findings),
            },
            "findings": [
                {
                    "severity": f.severity.value,
                    "category": f.category,
                    "feature": f.feature,
                    "detail": f.detail,
                    "suggestion": f.suggestion,
                }
                for f in self.gap.findings
            ],
        }
        return json.dumps(data, indent=2)

    def to_markdown(self) -> str:
        lines = ["# Gap Analysis Report", ""]
        lines.append("## Summary")
        lines.append("")
        lines.append("| Severity | Count |")
        lines.append("|----------|-------|")
        lines.append(f"| Supported / Mapped | {self.gap.n_supported()} |")
        lines.append(f"| Approximated | {self.gap.n_approximated()} |")
        lines.append(f"| Unsupported (non-blocking) | {self.gap.n_unsupported_nonblocking()} |")
        lines.append(f"| Unsupported (blocking) | {self.gap.n_unsupported_blocking()} |")
        lines.append(f"| Unavailable (source info missing) | {self.gap.n_unavailable()} |")
        lines.append(f"| **Total** | {len(self.gap.findings)} |")
        lines.append("")

        if self.gap.has_blocking():
            lines.append("**Blocking incompatibilities:** YES — conversion blocked")
        else:
            lines.append("**Blocking incompatibilities:** NONE")
        lines.append("")

        categories = [
            ("## 1. Supported / Mapped Features", Severity.SUPPORTED),
            ("## 2. Approximated / Fallback Mappings", Severity.APPROXIMATED),
            ("## 3. Unavailable Source Information", Severity.UNAVAILABLE),
            ("## 4. Unsupported Features (non-blocking)", Severity.UNSUPPORTED_NONBLOCK),
            ("## 5. Unsupported Features (blocking)", Severity.UNSUPPORTED_BLOCK),
        ]

        for header, sev in categories:
            items = [f for f in self.gap.findings if f.severity == sev]
            lines.append(header)
            lines.append("")
            if items:
                lines.append("| Category | Feature | Detail | Suggestion |")
                lines.append("|----------|---------|--------|------------|")
                for f in items:
                    lines.append(
                        f"| {f.category} | {f.feature} | {f.detail} | "
                        f"{f.suggestion if f.suggestion else '—'} |"
                    )
            else:
                lines.append("_None._")
            lines.append("")

        return "\n".join(lines)
