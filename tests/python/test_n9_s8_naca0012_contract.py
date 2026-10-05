"""N9-S8 NACA0012 validation-contract tests."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from import_tmr_naca0012 import infer_dimensions

def test_tmr_grid_filename_contract():
    assert infer_dimensions(Path("n0012_225-65.p2dfmt"))==(225,65)
    assert infer_dimensions(Path("n0012_449-129.p2dfmt"))==(449,129)
    assert infer_dimensions(Path("n0012_897-257.p2dfmt"))==(897,257)

def test_n9_s8_reference_document_exists():
    text=(ROOT/"docs/validation/N9-S8_NACA0012.md").read_text(encoding="utf-8")
    for required in ("NASA TMR","N9-S8","Re=1000","alpha=0","225x65","449x129","897x257","Total Validation","NOT QUALIFIED"):
        assert required in text
