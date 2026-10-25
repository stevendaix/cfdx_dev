#!/usr/bin/env python3
"""Qualify reading the SU2 NACA0012 Euler fixture, with explicit scope.

This fixture is not a CFDX aerodynamic validation.  Its compressible Euler
physics and JST discretization have no production path in the current
incompressible executable.  The faithful converter must keep reporting
unsupported numerical selections instead of silently replacing them.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "tests" / "data" / "su2"
sys.path.insert(0, str(ROOT / "src" / "cfdx" / "python"))

from cfdx.io.adapters.su2 import Su2Adapter  # noqa: E402
from cfdx.io.converter import convert  # noqa: E402
from cfdx.io.gap_analysis import GapAnalysisReport  # noqa: E402
from cfdx.io.interfaces import ConversionResult, SourceInfo  # noqa: E402


def qualify_reader(directory: Path) -> dict:
    """Retain separate adapter, canonical-conversion and physics verdicts."""
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(DATA_DIR / "mesh_NACA0012_inv.su2", directory / "mesh.su2")
    shutil.copy2(DATA_DIR / "inv_NACA0012_basic.cfg", directory / "config.cfg")
    result = ConversionResult(source=SourceInfo())
    read_ok = Su2Adapter().convert(directory, result)
    setup = result.setup
    checks = {
        "reader_completed": bool(read_ok),
        "no_blocking_reader_finding": not result.gap_report.has_blocking(),
        "points": result.mesh is not None and result.mesh["points"].shape == (5233, 3),
        "triangles": result.mesh is not None and len(result.mesh["cells"].get("triangle", [])) == 10216,
        "physics": setup is not None and setup.physics_model == "compressible_euler",
        "energy": setup is not None and setup.energy_model == "compressible",
        "boundaries": setup is not None and {bc.patch_name for bc in setup.boundary_conditions} == {"airfoil", "farfield"},
        "pressure": setup is not None and setup.initial_condition.pressure == 101325.0,
        "temperature": setup is not None and setup.initial_condition.temperature == 288.15,
    }
    # Keep the stricter converter result visible.  A reader PASS must not
    # be interpreted as acceptance of unsupported JST or compressible flow.
    canonical = convert(directory)
    boundary_by_name = {bc.patch_name: bc for bc in setup.boundary_conditions} if setup else {}
    airfoil = boundary_by_name.get("airfoil")
    farfield = boundary_by_name.get("farfield")
    velocity = setup.initial_condition if setup else None
    physical_metadata_checks = {
        "euler_wall_slip": airfoil is not None and airfoil.value_type.value == "wall_slip",
        "initial_velocity_vector_matches_speed": velocity is not None and math.isclose(
            math.sqrt(sum(value * value for value in velocity.velocity_vector)),
            velocity.velocity, rel_tol=1e-12, abs_tol=1e-12),
        "farfield_velocity_matches_initial_speed": velocity is not None and farfield is not None and math.isclose(
            math.sqrt(sum(value * value for value in farfield.velocity_vector)) if farfield.velocity_vector
            else farfield.velocity_magnitude,
            velocity.velocity, rel_tol=1e-12, abs_tol=1e-12),
    }
    return {
        "scope": "SU2 NACA0012 input-reader qualification",
        "reader_passed": all(checks.values()),
        "reader_checks": checks,
        "physical_metadata_checks": physical_metadata_checks,
        "boundary_initial_condition_equivalence_passed": all(physical_metadata_checks.values()),
        "reader_findings": json.loads(GapAnalysisReport(result.gap_report).to_json()),
        "canonical_conversion_passed": canonical.success,
        "canonical_conversion_findings": json.loads(GapAnalysisReport(canonical.gap_analysis).to_json()),
        "production_supported": False,
        "production_reason": "The CFDX production executable supports an incompressible carrier; this fixture declares compressible Euler at Mach 0.8.",
        "physical_validation": "NOT_RUN",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="directory for input copies and the reader report")
    parser.add_argument("--reader-only", action="store_true",
                        help="request only format-reader checks; permit success without a flow solve")
    args = parser.parse_args(argv)
    directory = args.output_dir or Path(tempfile.mkdtemp(prefix="cfdx_naca0012_reader_"))
    report = qualify_reader(directory)
    report_path = directory / "reader_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("=== SU2 NACA0012 input-reader qualification ===")
    print("Reader:", "PASS" if report["reader_passed"] else "FAIL")
    print("Boundary/initial-condition equivalence:",
          "PASS" if report["boundary_initial_condition_equivalence_passed"] else "FAIL")
    for check, passed in report["physical_metadata_checks"].items():
        if not passed:
            print(f"  MISMATCH: {check}")
    print("Canonical conversion:", "PASS" if report["canonical_conversion_passed"] else "BLOCKED")
    for finding in report["canonical_conversion_findings"]["findings"]:
        if finding["severity"] == "unsupported_blocking":
            print(f"  BLOCKING: {finding['feature']} - {finding['detail']}")
    print("Production: UNSUPPORTED for this compressible Euler fixture")
    print("Physical validation: NOT_RUN; no Cd, Cl or Cp prediction is asserted")
    print(f"Reader evidence: {report_path}")
    if not report["reader_passed"]:
        return 1
    if args.reader_only:
        return 0
    # Preserve a failing default verdict for legacy callers that expected
    # this command to validate the full aerodynamic run.
    print("Full aerodynamic validation is unavailable. Use --reader-only for the explicitly limited reader scope.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
