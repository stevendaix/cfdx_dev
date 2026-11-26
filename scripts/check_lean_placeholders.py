"""Reject explicit proof placeholders and local axioms in CFDX Lean sources."""

from __future__ import annotations

import re
import sys
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[1] / "formal" / "lean" / "CFDX"
FORBIDDEN = {
    "sorry": re.compile(r"\bsorry\b"),
    "admit": re.compile(r"\badmit\b"),
    # Lean permits declaration modifiers before axiom (notably private).
    # Check the declaration command after supported command modifiers.
    "axiom declaration": re.compile(
        r"(?m)^\s*(?:(?:private|protected|noncomputable|unsafe|partial|scoped|local|global)\s+)*"
        r"axiom\s+[A-Za-z_][A-Za-z0-9_'.]*"
    ),
}


def find_violations(source: str, display_path: str) -> list[str]:
    """Return policy violations found in one Lean source string."""
    violations: list[str] = []
    for label, pattern in FORBIDDEN.items():
        for match in pattern.finditer(source):
            line = source.count("\n", 0, match.start()) + 1
            violations.append(f"{display_path}:{line}: {label}")
    return violations


def main() -> int:
    if not SOURCE_ROOT.is_dir():
        print(f"Lean source directory not found: {SOURCE_ROOT}", file=sys.stderr)
        return 1

    violations: list[str] = []
    for path in sorted(SOURCE_ROOT.rglob("*.lean")):
        source = path.read_text(encoding="utf-8")
        violations.extend(
            find_violations(source, str(path.relative_to(SOURCE_ROOT.parent.parent)))
        )

    if violations:
        print("Forbidden unproved declarations/placeholders found:")
        print("\n".join(violations))
        return 1

    print("Lean source placeholder check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
