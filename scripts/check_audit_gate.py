#!/usr/bin/env python3
"""Audit gate for CFDX (constitution technique, Rules of Gold 1 and 2).

Modes
-----
1. --validate-fiches
   Structural validation of every docs/audits/*.md fiche methode:
   required sections, evidence header, at least one DOI/ISBN reference.
2. --changed-files <file-list> [--pr-body <file>]
   Rule of Gold 1: any change under src/ (numerical core) must be accompanied
   by a reference to an approved fiche in docs/audits/ (either a changed
   docs/audits/*.md file in the same PR, or an explicit docs/audits/... path
   mentioned in the PR body).

Exit code 1 on any violation (CI gate). Pure stdlib, Python >= 3.9.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REQUIRED_SECTIONS = [
    "## 1. R\u00e9f\u00e9rences sources",
    "## 2. Formulation continue",
    "## 3. Discr\u00e9tisation Volume Fini",
    "## 4. Traitement des fronti\u00e8res",
    "## 5. Pi\u00e8ges connus",
    "## 6. Ordre de pr\u00e9cision",
    "## 7. Statut d'impl\u00e9mentation",
]

REF_PATTERN = re.compile(r"(10\.\d{4,9}/\S+|ISBN[\s:]*[0-9\-\u2013 ]{9,17})", re.IGNORECASE)
NUMERIC_PREFIXES = ("src/cfdx/core/numerics/", "src/cfdx/core/linalg/", "src/cfdx/physics/")


def validate_fiches(root: Path) -> list[str]:
    errors: list[str] = []
    audit_dir = root / "docs" / "audits"
    fiches = sorted(p for p in audit_dir.glob("*.md") if not p.name.startswith("README") and not p.name.startswith("_"))
    if not fiches:
        return ["docs/audits contains no fiche methode"]
    for fiche in fiches:
        text = fiche.read_text(encoding="utf-8")
        for section in REQUIRED_SECTIONS:
            if section not in text:
                errors.append(f"{fiche.name}: missing section '{section}'")
        if not REF_PATTERN.search(text):
            errors.append(f"{fiche.name}: no DOI or ISBN reference found (Rule of Gold 1)")
        if "Statut\u00a0:" not in text and "Statut:" not in text:
            errors.append(f"{fiche.name}: missing evidence header 'Statut'")
    return errors


def check_audit_gate(root: Path, changed: list[str], body: str) -> list[str]:
    errors: list[str] = []
    numeric = [c for c in changed if c.startswith(NUMERIC_PREFIXES) and c.endswith((".h", ".cpp", ".cu", ".hpp"))]
    if not numeric:
        return errors
    audit_files = [c for c in changed if c.startswith("docs/audits/")]
    body_refs = re.findall(r"docs/audits/[A-Za-z0-9_.\u2013-]+\.md", body)
    if not audit_files and not body_refs:
        errors.append(
            "Rule of Gold 1 violation: this PR modifies the numerical core ("
            + ", ".join(numeric[:3])
            + (" ..." if len(numeric) > 3 else "")
            + ") but references no fiche methode in docs/audits/. "
              "Add or update a fiche, or link an approved one in the PR description."
        )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="CFDX audit gate")
    parser.add_argument("--root", default=".", help="repository root")
    parser.add_argument("--validate-fiches", action="store_true")
    parser.add_argument("--changed-files", help="file containing one changed path per line")
    parser.add_argument("--pr-body", help="file containing the PR description")
    args = parser.parse_args()

    root = Path(args.root)
    errors: list[str] = []

    if args.validate_fiches:
        errors.extend(validate_fiches(root))

    if args.changed_files:
        changed = [line.strip() for line in Path(args.changed_files).read_text(encoding="utf-8").splitlines() if line.strip()]
        body = Path(args.pr_body).read_text(encoding="utf-8") if args.pr_body else ""
        errors.extend(check_audit_gate(root, changed, body))

    if errors:
        print("AUDIT GATE FAILED")
        for e in errors:
            print("  - " + e)
        return 1
    print("AUDIT GATE PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
