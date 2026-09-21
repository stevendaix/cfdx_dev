#!/usr/bin/env python3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cfdx import Case


def main() -> int:
    case = (
        Case("api-test")
        .enable("incompressible")
        .set_numerics(gradient="gauss_linear", laplacian="corrected")
        .set_boundary("inlet", velocity=[1.0, 0.0, 0.0])
    )
    data = case.as_dict()
    assert data["name"] == "api-test"
    assert data["physics"]["incompressible"]["enabled"] is True
    assert data["numerics"]["gradient"] == "gauss_linear"
    assert data["boundaries"]["inlet"]["velocity"][0] == 1.0
    print("CFDX_PYTHON_API_TEST: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
