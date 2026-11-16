"""N9-S8 NACA0012 validation-contract tests; no pytest dependency required."""

import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT / "scripts"))
sys.path.append(str(ROOT / "scripts"))
from import_tmr_naca0012 import infer_dimensions


def main():
    assert infer_dimensions(Path("n0012_225-65.p2dfmt")) == (225, 65)
    assert infer_dimensions(Path("n0012_449-129.p2dfmt")) == (449, 129)
    assert infer_dimensions(Path("n0012_897-257.p2dfmt")) == (897, 257)
    text = (ROOT / "docs/validation/14-benchmarks/06-naca0012/chapter.py").read_text(
        encoding="utf-8"
    )
    # Jupytext markdown-only chapters remain importable Python source.
    # The prior escaped LaTeX rho/nu split had created uncompilable code.
    ast.parse(text, filename="docs/validation/14-benchmarks/06-naca0012/chapter.py")
    for required in (
        "NASA/TMR",
        "N9-S8",
        "Re=1000",
        "alpha=0",
        "225 x 65",
        "449 x 129",
        "897 x 257",
        "Total Validation",
        "NOT QUALIFIED",
    ):
        assert required in text, required
    print("N9-S8-NACA0012-CONTRACT: PASS")


if __name__ == "__main__":
    main()
