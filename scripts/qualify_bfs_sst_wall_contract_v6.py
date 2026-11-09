"""Fail-closed BFS pressure/SST/wall-function comparability preflight.

Never treat a converged pressure *substep* as a converged RANS simulation.
OpenFOAM's reference wall functions are not interchangeable with CFDX's
experimental fixed-k/wall-omega-factor boundary treatment.  Native CFDX
wall CSV and the independent three-grid physical oracles are still required.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

REFERENCE = {
    "source": "OpenFOAM 2606 backwardFacingStep2D / Driver-Seegmiller",
    "url": "https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/",
    "h_m": 0.0127,
    "u_inlet_m_s": 44.2,
    "re_h": 36000,
    "reattachment_x_over_h": 6.26,
    "reattachment_uncertainty_absolute": 0.10,
    "openfoam_wall_nut": "nutUSpaldingWallFunction",
    "openfoam_wall_k": "kqRWallFunction",
    "openfoam_wall_omega": "omegaWallFunction",
}


class BFSWallContractError(ValueError):
    pass


def load_case(hdf5_path: Path) -> dict[str, Any]:
    try:
        import h5py  # type: ignore[import-not-found]

        with h5py.File(hdf5_path) as stream:
            raw = stream.attrs.get("case_setup_json")
        if not isinstance(raw, (str, bytes)):
            raise BFSWallContractError("missing case_setup_json in native HDF5")
        case = json.loads(raw)
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        raise BFSWallContractError(f"cannot load case metadata: {exc}") from exc
    if not isinstance(case, dict):
        raise BFSWallContractError("invalid case object")
    return case


def inspect(case: dict[str, Any], log: str = "", convergence: Any = None) -> dict:
    if not isinstance(case, dict):
        raise BFSWallContractError("invalid case JSON")
    boundaries = case.get("boundary_conditions")
    numerics = case.get("numerics")
    if not isinstance(boundaries, list) or not isinstance(numerics, dict):
        raise BFSWallContractError("missing native boundary/numerical contracts")
    raw = numerics.get("raw_settings", {})
    if not isinstance(raw, dict):
        raise BFSWallContractError("invalid raw numerical settings")
    walls = [b for b in boundaries if isinstance(b, dict) and b.get("type") == "wall"]
    if not walls:
        raise BFSWallContractError("no explicit solid wall patches")
    sst_methods = {
        b.get("patch_name", "<unnamed>"): {
            s.get("field"): s.get("type")
            for s in b.get("sst_scalar_conditions", [])
            if isinstance(s, dict)
        }
        for b in walls
    }
    config = {
        "wall_omega_factor": raw.get("sst_wall_omega_factor", "6"),
        "sst_molecular_wall_diffusion": raw.get(
            "sst_molecular_wall_diffusion", "false"
        ),
        "momentum_scheme": numerics.get("momentum_scheme", "unknown"),
        "turbulence_model": case.get("turbulence_model", "unknown"),
        "sst_boundary_methods": sst_methods,
    }
    # The names in a native CFDX case describe its own experimental BCs,
    # NOT an implemented OpenFOAM wall-function contract. No string alias
    # or nominal omega factor constitutes an equivalence proof.
    equivalent_wall_functions = False
    reasons = [
        ("CFDX does not provide verified equivalents of OpenFOAM "
         "nutUSpaldingWallFunction/kqRWallFunction/omegaWallFunction for this case")
    ]
    if raw.get("sst_molecular_wall_diffusion", "false") != "true":
        reasons.append("molecular near-wall diffusion path is legacy/disabled")
    if numerics.get("momentum_scheme") == "first_order":
        reasons.append(
            "convection is first order rather than independently qualified"
        )
    pressure_records = []
    for line in log.splitlines():
        if not line.startswith("PRESSURE_HILO_PHYSICAL_FLUX "):
            continue
        fields = dict(re.findall(r"([a-z_][a-z_0-9]*)=([^\s]+)", line))
        try:
            physical = float(fields["exact_relative"])
            rounded = float(fields["stored_fp64_relative"])
        except (KeyError, ValueError) as exc:
            raise BFSWallContractError("malformed physical pressure record") from exc
        if not math.isfinite(physical) or not math.isfinite(rounded):
            raise BFSWallContractError("nonfinite physical pressure record")
        pressure_records.append(
            {
                "exact_relative": physical,
                "stored_fp64_relative": rounded,
                "substep_pass": physical <= 1e-11 and rounded <= 1e-11,
            }
        )
    if not isinstance(convergence, (dict, type(None))):
        raise BFSWallContractError("invalid convergence certificate")
    nonlinear_converged = (
        convergence is not None and convergence.get("converged") is True
    )
    pressure_substeps_pass = bool(pressure_records) and all(
        r["substep_pass"] for r in pressure_records
    )
    # Independent quantitative evidence must exist; even if the solver
    # reports converged, this audit cannot silently certify wall comparability.
    return {
        "case_id": "TH-03-BFS",
        "qualification_status": "BLOCKED_MODEL_AND_PHYSICAL_EVIDENCE",
        "reference": REFERENCE,
        "cfdx_wall_contract": config,
        "pressure_substep_records": pressure_records,
        "pressure_substeps_pass_original_1e_11": pressure_substeps_pass,
        "solver_nonlinear_converged": nonlinear_converged,
        "openfoam_wall_functions_equivalent": equivalent_wall_functions,
        "physical_qualification": False,
        "reattachment_qualified": False,
        "wall_yplus_qualified": False,
        "wall_cf_qualified": False,
        "velocity_profiles_qualified": False,
        "reasons": reasons
        + [
            "need accepted converged momentum, continuity, SST k/omega budgets",
            "need actual native wall-face y+/Cf and independent velocity profiles",
            "need wall-model equivalence or a separately validated resolved-wall comparison",
            "need mesh sensitivity and physical/experimental references",
        ],
        "available_next_oracles": [
            "scripts/audit_th03_wall_resolution.py",
            "scripts/extract_th03_wall_diagnostics.py",
            "scripts/qualify_thermohydraulic_benchmarks.py",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-hdf5", type=Path, required=True)
    parser.add_argument("--run-log", type=Path)
    parser.add_argument("--convergence-json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        case = load_case(args.case_hdf5)
        log = args.run_log.read_text(encoding="utf-8") if args.run_log else ""
        convergence = (
            json.loads(args.convergence_json.read_text(encoding="utf-8"))
            if args.convergence_json
            else None
        )
        result = inspect(case, log, convergence)
    except (BFSWallContractError, OSError, ValueError) as exc:
        parser.exit(2, f"BFS wall/SST gate ERROR: {exc}\n")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(
        "BFS wall/SST: model non-equivalence and missing physical qualification; BLOCKED"
    )
    # This is a certification gate, not a preflight which may exit green.
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
