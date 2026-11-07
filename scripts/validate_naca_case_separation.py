"""Fail-closed separation of CFDX laminar NACA and NASA/TMR SST campaigns.

This checks case declarations and *optional* TMR mesh provenance. It does NOT
qualify a turbulent flow or infer an SST solution from imported mesh metadata.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def check_registry(data: dict) -> tuple[dict, dict]:
    cases = data.get("cases")
    if not isinstance(cases, list):
        raise TypeError("case registry must be a list")
    by_id = {c.get("id"): c for c in cases if isinstance(c, dict)}
    if len(by_id) != len(cases):
        raise ValueError("duplicate/invalid NACA qualification case IDs")
    laminar = by_id.get("EXT-NACA0012")
    sst = by_id.get("TUR-NACA0012")
    if laminar is None or sst is None:
        raise ValueError("both EXT-NACA0012 and TUR-NACA0012 must exist")
    if laminar is sst or laminar.get("domain") == sst.get("domain"):
        raise ValueError("laminar and turbulent NACA scopes are conflated")
    if "laminar" not in str(laminar.get("solver_settings", "")).lower():
        raise ValueError("EXT-NACA0012 must retain its laminar subcase")
    if "sst" not in str(sst.get("numerical_schemes", "")).lower():
        raise ValueError("TUR-NACA0012 must explicitly name SST")
    if "re=6e6" not in str(sst.get("parameters", "")).lower().replace(" ", ""):
        raise ValueError("TUR-NACA0012 must retain its Re=6e6 contract")
    # This tool does not read a physical solution or independent reference.
    # No PASS can be certified using only the registry or mesh metadata.
    if sst.get("status") == "PASS":
        raise ValueError("TUR-NACA0012 PASS requires independent physical evidence, not this metadata gate")
    return laminar, sst


def check_tmr_mesh(path: Path) -> str:
    import h5py

    with h5py.File(path, "r") as h:
        if h.attrs.get("validation_source") != "NASA TMR NACA0012":
            raise ValueError("SST case requires explicitly sourced NASA/TMR NACA mesh")
        if h.attrs.get("validation_status") != "imported-provenance-only":
            raise ValueError("TMR mesh must remain import-only, not a fabricated CFD result")
        if "turbmodels/naca0012_grids.html" not in str(
            h.attrs.get("validation_source_url", "")
        ):
            raise ValueError("unexpected NACA mesh family source")
        fingerprint = str(h.attrs.get("validation_source_sha256", ""))
        if len(fingerprint) != 64 or any(ch not in "0123456789abcdef" for ch in fingerprint):
            raise ValueError("missing/invalid TMR source sha256")
        patch_types = {}
        for entry in str(h.attrs.get("boundary_patches", "")).split(";"):
            if entry:
                parts = entry.split(":")
                if len(parts) != 4:
                    raise ValueError("malformed patch metadata")
                patch_types[parts[0]] = int(parts[-1])
        for name, kind in {"airfoil": 0, "front": 6, "back": 6}.items():
            if patch_types.get(name) != kind:
                raise ValueError(f"TMR {name} patch contract mismatch")
        if "points" not in h or "owner" not in h or "neighbour" not in h:
            raise ValueError("incomplete TMR mesh interchange format")
        n_cells = len(h["cell_offsets"]) - 1
        if n_cells <= 0:
            raise ValueError("empty TMR mesh")
        return str(h.attrs.get("validation_grid", "unknown"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=Path(
        "docs/validation/CFDX_QUALIFICATION_REGISTRY.json"))
    parser.add_argument("--tmr-mesh", type=Path, help="optional imported NASA/TMR HDF5")
    args = parser.parse_args()
    try:
        _, sst = check_registry(json.loads(args.registry.read_text(encoding="utf-8")))
        grid = check_tmr_mesh(args.tmr_mesh) if args.tmr_mesh else "not-provided"
    except (OSError, ValueError, KeyError) as exc:
        print(f"NACA_CASE_SEPARATION: FAIL {exc}")
        return 2
    print(f"NACA_CASE_SEPARATION: PASS laminar=EXT-NACA0012 "
          f"turbulent=TUR-NACA0012 sst_status={sst['status']} "
          f"tmr_mesh={grid} scope=METADATA_ONLY_NOT_PHYSICAL_VALIDATION")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
