"""Regression tests for the Lean source proof-policy guard."""

from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SCRIPT_PATH = (
    Path(__file__).resolve().parents[2] / "scripts" / "check_lean_placeholders.py"
)
SPEC = spec_from_file_location("check_lean_placeholders", SCRIPT_PATH)
assert SPEC is not None
assert SPEC.loader is not None
GUARD = module_from_spec(SPEC)
SPEC.loader.exec_module(GUARD)


@pytest.mark.parametrize(
    "declaration",
    [
        "axiom bogus : False",
        "private axiom bogus : False",
        "protected axiom bogus : False",
        "noncomputable private axiom bogus : False",
        "unsafe axiom bogus : False",
        "scoped private axiom bogus : False",
        "local axiom bogus : False",
        "global axiom bogus : False",
    ],
)
def test_axiom_declarations_with_modifiers_are_rejected(declaration: str) -> None:
    violations = GUARD.find_violations(declaration, "fixture.lean")
    assert len(violations) == 1
    assert "axiom declaration" in violations[0]


@pytest.mark.parametrize("placeholder", ["sorry", "admit"])
def test_unproved_placeholders_are_rejected(placeholder: str) -> None:
    violations = GUARD.find_violations(f"theorem bogus : True := {placeholder}", "fixture.lean")
    assert len(violations) == 1
    assert placeholder in violations[0]


def test_proved_theorem_is_accepted() -> None:
    assert GUARD.find_violations("theorem good : True := by trivial", "fixture.lean") == []
