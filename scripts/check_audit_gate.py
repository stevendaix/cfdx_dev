"""Audit gate for CFDX (constitution technique, Rules of Gold 1 and 2).

Modes
-----
1. --validate-fiches
   Structural validation of every docs/audits/*.md fiche methode:
   numbered sections 1..7, evidence header, at least one DOI/ISBN reference.
   Numerical-method fiches must additionally carry the exact math-section
   titles (2: formulation continue, 3: discretisation Volume Fini, 4:
   traitement des frontieres, 5: pieges connus). Non-method fiches
   (oracles/benchmark catalogue, reference-code policy) are exempt from the
   math-title check but keep the numbered structure, Statut header and
   DOI/ISBN requirements.
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

# Some fiche headers use a typographic non-breaking space (U+00A0) before
# the colon (evidence header and ISBN lines). Normalise it to a plain
# space before every structural check so both spellings are accepted.
NBSP = chr(0xA0)

# Structural invariant shared by every fiche: sections numbered 1..7.
REQUIRED_NUMBERED_SECTIONS = tuple(f"## {i}." for i in range(1, 8))

# Exact math-contract titles required on numerical-method fiches only.
MATH_SECTION_TITLES = (
    "## 2. Formulation continue",
    "## 3. Discrétisation Volume Fini",
    "## 4. Traitement des frontières",
    "## 5. Pièges connus",
)
# Section 7 title accepts both apostrophe variants (ASCII and typographic).
SECTION7_VARIANTS = (
    "## 7. Statut d'implémentation",
    "## 7. Statut d’implémentation",
)
# Fiches that are not numerical-method contracts (oracles catalogue,
# reference-code policy): exempt from math-title checks, not from structure.
NON_METHOD_FICHES = {
    "07-benchmarks-references.md",
    "08-codes-de-reference.md",
}

# Consolidated audit dossiers and active tracking registers are reports, not
# numerical-method fiches. They retain historical sections and checklists and
# therefore do not follow the strict 1..7 method-fiche template.
NON_FICHE_REPORTS = {
    "audit_v2.md",
    "turbulence.md",
    # 2026-10-25 consolidated correction/integration dossiers (master 026cb499):
    # review reports without the method-fiche template, same category as
    # audit_v2.md. Excluded from structural fiche validation only.
    "CFDX_CODE_CORRECTIONS_COMPLET_2026-10-25.md",
    "CFDX_CODE_ET_VALIDATION_6_VERROUX_2026-10-25.md",
    "CFDX_etat_integration_et_reste_a_faire_2026-10-25.md",
}

# Reference detection (DOI prefix or ISBN). Every character class avoids
# backslash escapes ([.] for the dot, [0-9] for digits, etc.) so the
# pattern cannot be damaged again by transport layers that strip raw
# escape characters; matching behaviour is identical to the escaped form.
REF_PATTERN = re.compile(
    r"(10[.][0-9]{4,9}/[0-9A-Za-z._()/-]+|ISBN[s: ]*[-0-9][0-9 -]{8,16})",
    re.IGNORECASE,
)
NUMERIC_PREFIXES = (
    "src/cfdx/core/numerics/",
    "src/cfdx/core/linalg/",
    "src/cfdx/physics/",
)


def validate_fiches(root: Path) -> list[str]:
    errors: list[str] = []
    audit_dir = root / "docs" / "audits"
    fiches = sorted(
        p
        for p in audit_dir.glob("*.md")
        if not p.name.startswith("README") and not p.name.startswith("_")
    )
    if not fiches:
        return ["docs/audits contains no fiche methode"]
    for fiche in fiches:
        if fiche.name in NON_FICHE_REPORTS:
            continue
        text = fiche.read_text(encoding="utf-8").replace(NBSP, " ")
        for section in REQUIRED_NUMBERED_SECTIONS:
            if section not in text:
                errors.append(f"{fiche.name}: missing numbered section '{section}'")
        if fiche.name not in NON_METHOD_FICHES:
            for title in MATH_SECTION_TITLES:
                if title not in text:
                    errors.append(f"{fiche.name}: missing math section '{title}'")
        if not any(variant in text for variant in SECTION7_VARIANTS):
            errors.append(f"{fiche.name}: missing section 7 'Statut d'implémentation'")
        if not REF_PATTERN.search(text):
            errors.append(
                f"{fiche.name}: no DOI or ISBN reference found (Rule of Gold 1)"
            )
        if "Statut :" not in text and "Statut:" not in text:
            errors.append(f"{fiche.name}: missing evidence header 'Statut'")
    return errors


def check_audit_gate(root: Path, changed: list[str], body: str) -> list[str]:
    errors: list[str] = []
    numeric = [
        c
        for c in changed
        if c.startswith(NUMERIC_PREFIXES) and c.endswith((".h", ".cpp", ".cu", ".hpp"))
    ]
    if not numeric:
        return errors
    audit_files = [c for c in changed if c.startswith("docs/audits/")]
    body_refs = re.findall(r"docs/audits/[A-Za-z0-9_.–-]+.md", body)
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
    parser.add_argument(
        "--changed-files", help="file containing one changed path per line"
    )
    parser.add_argument("--pr-body", help="file containing the PR description")
    args = parser.parse_args()

    root = Path(args.root)
    errors: list[str] = []

    if args.validate_fiches:
        errors.extend(validate_fiches(root))

    if args.changed_files:
        raw_changed = Path(args.changed_files).read_text(encoding="utf-8")
        changed = [line.strip() for line in raw_changed.splitlines() if line.strip()]
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
