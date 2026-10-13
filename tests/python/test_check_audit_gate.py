"""Contract tests for scripts/check_audit_gate.py (Rules of Gold 1 and 2).

The audit gate is a blocking CI control. These tests verify that it really
rejects the forbidden cases -- malformed fiches, numerical-core changes
without an approved fiche reference, and invalid PR bodies -- instead of
only ever checking the happy path (the negative cases requested by the
2026-10-12 audit of the 10 last merged PRs).

All synthetic fiche content is built from the script's own constants so
these tests cannot drift from the exact (accented) section titles.
"""

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "check_audit_gate.py"
SPEC = importlib.util.spec_from_file_location("cfdx_check_audit_gate", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)

FICHE_NAME = "01-gate-fixture.md"
NUMERIC_CORE = "src/cfdx/physics/wall_distance.h"
DOI_LINE = "DOI: 10.1016/0045-7825(74)90029-2"


def valid_fiche_text():
    lines = [
        "# FICHE METHODE (synthetic audit-gate test fixture)",
        "",
        "Statut : APPROUVEE",
        "",
        "## 1. References",
        *GATE.MATH_SECTION_TITLES,
        "## 6. Precision and validation",
        GATE.SECTION7_VARIANTS[0],
        "",
        DOI_LINE,
        "",
    ]
    return "\n".join(lines)


def write_fiche(root, name=FICHE_NAME, text=None):
    audits = root / "docs" / "audits"
    audits.mkdir(parents=True, exist_ok=True)
    if text is None:
        text = valid_fiche_text()
    (audits / name).write_text(text, encoding="utf-8")
    return audits / name


def run_gate_cli(root, changed=None, body=None, validate_fiches=False):
    cmd = [sys.executable, str(SCRIPT), "--root", str(root)]
    if validate_fiches:
        cmd.append("--validate-fiches")
    if changed is not None:
        cmd += ["--changed-files", str(changed)]
    if body is not None:
        cmd += ["--pr-body", str(body)]
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


class ValidateFichesContractTests(unittest.TestCase):
    """Rule of Gold 1: structural validation of the fiches methode."""

    def test_repository_fiches_are_structurally_valid(self):
        self.assertEqual(GATE.validate_fiches(ROOT), [])

    def test_empty_audit_directory_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            audits = Path(tmp) / "docs" / "audits"
            audits.mkdir(parents=True)
            self.assertEqual(
                GATE.validate_fiches(Path(tmp)),
                ["docs/audits contains no fiche methode"],
            )

    def test_valid_fiche_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fiche(Path(tmp))
            self.assertEqual(GATE.validate_fiches(Path(tmp)), [])

    def test_missing_numbered_section_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = valid_fiche_text().replace("## 6. Precision and validation\n", "")
            write_fiche(Path(tmp), text=text)
            errors = GATE.validate_fiches(Path(tmp))
            self.assertTrue(
                any("missing numbered section '## 6.'" in e for e in errors)
            )

    def test_missing_math_title_is_rejected_on_method_fiche(self):
        with tempfile.TemporaryDirectory() as tmp:
            title = GATE.MATH_SECTION_TITLES[1]
            text = valid_fiche_text().replace(title + "\n", "")
            write_fiche(Path(tmp), text=text)
            errors = GATE.validate_fiches(Path(tmp))
            self.assertTrue(any("missing math section" in e for e in errors))
            self.assertTrue(
                any("missing numbered section '## 3.'" in e for e in errors)
            )

    def test_non_method_fiche_is_exempt_from_math_titles_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            title = GATE.MATH_SECTION_TITLES[1]
            # A non-method fiche is exempt from the mandatory math
            # titles, but not from the numbered structure: swap the math
            # title for a plain "## 3." section instead of deleting
            # the line, otherwise the structural check rejects the fixture.
            text = valid_fiche_text().replace(title, "## 3. Benchmarks and references")
            write_fiche(Path(tmp), name="07-benchmarks-references.md", text=text)
            self.assertEqual(GATE.validate_fiches(Path(tmp)), [])

    def test_fiche_without_reference_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = valid_fiche_text().replace(DOI_LINE + "\n", "")
            write_fiche(Path(tmp), text=text)
            errors = GATE.validate_fiches(Path(tmp))
            self.assertTrue(any("no DOI or ISBN reference" in e for e in errors))

    def test_fiche_without_statut_header_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = valid_fiche_text().replace("Statut : APPROUVEE\n", "")
            write_fiche(Path(tmp), text=text)
            errors = GATE.validate_fiches(Path(tmp))
            self.assertTrue(any("missing evidence header" in e for e in errors))

    def test_missing_section7_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = valid_fiche_text().replace(GATE.SECTION7_VARIANTS[0] + "\n", "")
            write_fiche(Path(tmp), text=text)
            errors = GATE.validate_fiches(Path(tmp))
            self.assertTrue(any("missing section 7" in e for e in errors))
            self.assertTrue(
                any("missing numbered section '## 7.'" in e for e in errors)
            )

    def test_typographic_section7_variant_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = valid_fiche_text().replace(
                GATE.SECTION7_VARIANTS[0], GATE.SECTION7_VARIANTS[1]
            )
            write_fiche(Path(tmp), text=text)
            self.assertEqual(GATE.validate_fiches(Path(tmp)), [])


class RuleOfGoldOneContractTests(unittest.TestCase):
    """Rule of Gold 1: a numerical-core change must cite an approved fiche."""

    def test_numeric_change_without_fiche_reference_is_rejected(self):
        errors = GATE.check_audit_gate(Path("/nonexistent"), [NUMERIC_CORE], "")
        self.assertEqual(len(errors), 1)
        self.assertIn("Rule of Gold 1 violation", errors[0])
        self.assertIn(NUMERIC_CORE, errors[0])

    def test_numeric_change_with_changed_fiche_is_accepted(self):
        changed = [NUMERIC_CORE, "docs/audits/11-rayonnement.md"]
        errors = GATE.check_audit_gate(Path("/nonexistent"), changed, "")
        self.assertEqual(errors, [])

    def test_numeric_change_with_fiche_reference_in_body_is_accepted(self):
        body = "Audit reference: docs/audits/10-thermique-cht.md (approved)."
        errors = GATE.check_audit_gate(Path("/nonexistent"), [NUMERIC_CORE], body)
        self.assertEqual(errors, [])

    def test_invalid_fiche_reference_in_body_is_rejected(self):
        # A reference that is not a docs/audits/*.md path does not satisfy
        # Rule of Gold 1.
        body = "Audited against docs/audits/09-notes.txt (see issue #480)."
        errors = GATE.check_audit_gate(Path("/nonexistent"), [NUMERIC_CORE], body)
        self.assertEqual(len(errors), 1)

    def test_non_numeric_change_does_not_require_a_fiche(self):
        changed = [
            "src/cfdx/io/vtu/vtu_writer.cpp",
            "tests/unit/test_field.cpp",
            "README.md",
        ]
        errors = GATE.check_audit_gate(Path("/nonexistent"), changed, "")
        self.assertEqual(errors, [])

    def test_all_numeric_source_extensions_are_gated(self):
        for extension in (".h", ".cpp", ".cu", ".hpp"):
            path = "src/cfdx/physics/implementation" + extension
            errors = GATE.check_audit_gate(Path("/nonexistent"), [path], "")
            self.assertEqual(len(errors), 1)

    def test_non_source_files_under_numeric_prefixes_are_not_gated(self):
        path = "src/cfdx/physics/notes.txt"
        errors = GATE.check_audit_gate(Path("/nonexistent"), [path], "")
        self.assertEqual(errors, [])


class CliExitCodeContractTests(unittest.TestCase):
    """End-to-end: the gate must exit 1 on violations, 0 on a clean PR."""

    def test_cli_fails_when_numerical_core_has_no_fiche(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            changed = root / "changed.txt"
            changed.write_text("src/cfdx/core/numerics/flux.h\n", encoding="utf-8")
            proc = run_gate_cli(root, changed=changed)
            self.assertEqual(proc.returncode, 1)
            self.assertIn("AUDIT GATE FAILED", proc.stdout)
            self.assertIn("Rule of Gold 1 violation", proc.stdout)

    def test_cli_passes_on_a_well_formed_pr(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_fiche(root)
            fiche_ref = "docs/audits/" + FICHE_NAME
            changed_lines = ["src/cfdx/physics/energy_solver.h", fiche_ref]
            changed = root / "changed.txt"
            changed.write_text("\n".join(changed_lines) + "\n", encoding="utf-8")
            body = root / "body.txt"
            body.write_text("Fiche: " + fiche_ref, encoding="utf-8")
            proc = run_gate_cli(root, changed=changed, body=body, validate_fiches=True)
            self.assertEqual(proc.returncode, 0)
            self.assertIn("AUDIT GATE PASSED", proc.stdout)

    def test_cli_fails_when_a_fiche_is_structurally_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            text = valid_fiche_text().replace(DOI_LINE + "\n", "")
            write_fiche(root, text=text)
            proc = run_gate_cli(root, validate_fiches=True)
            self.assertEqual(proc.returncode, 1)
            self.assertIn("no DOI or ISBN reference", proc.stdout)


if __name__ == "__main__":
    unittest.main()
