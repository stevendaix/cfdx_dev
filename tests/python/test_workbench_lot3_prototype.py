"""Static and JavaScript checks for the optional real-data-only UX concept."""

from html.parser import HTMLParser
from pathlib import Path
import shutil
import subprocess

import pytest

PROTOTYPE = Path(__file__).resolve().parents[2] / "prototypes/gui-nextgen/modern_lot3.html"


class _ElementIds(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.ids: set[str] = set()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if attributes.get("id"):
            assert attributes["id"] not in self.ids, "duplicate HTML id"
            self.ids.add(attributes["id"])


def test_prototype_is_default_black_without_fabricated_solver_evidence() -> None:
    html = PROTOTYPE.read_text(encoding="utf-8")
    assert '<html lang="fr" data-theme="black">' in html
    assert "--bg:#000000" in html
    assert 'value="black"' in html
    assert "no current-HEAD" not in html  # Nothing masquerades as qualified
    assert "Aucune donnée chargée" in html
    assert "Non vérifié" in html
    parser = _ElementIds()
    parser.feed(html)
    assert {"themeSelect", "commandPalette", "commandQuery", "commandResults"} <= parser.ids
    # monitorFile is only created on the Run page; the prototype is still a shell.
    assert "const themes=['black','graphite','light']" in html
    assert "function searchCommands" in html
    assert "prompt(" not in html
    assert "addEventListener('keydown'" in html
    assert "document.createElement('td')" in html  # imported registry is textContent-only


@pytest.mark.skipif(shutil.which("node") is None, reason="Node.js optional")
def test_prototype_javascript_parses_without_a_browser() -> None:
    html = PROTOTYPE.read_text(encoding="utf-8")
    script = html.split("<script>", 1)[1].split("</script>", 1)[0]
    result = subprocess.run(
        ["node", "--check", "--input-type=commonjs"],
        input=script, text=True, capture_output=True, timeout=10, check=False,
    )
    assert result.returncode == 0, result.stderr
