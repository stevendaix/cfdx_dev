#!/usr/bin/env python3
"""CFDX I/O converter CLI.

Usage:
    cfdx-io convert <case_path> [--output DIR] [--solver NAME]
    cfdx-io list-adapters
"""
from __future__ import annotations

import sys

from cfdx.io.pipeline import convert, detect_adapter, get_available_adapters


def main(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]

    if not args:
        print(__doc__)
        return 1

    cmd = args[0]

    if cmd == "list-adapters":
        adapters = get_available_adapters()
        print("Available solver adapters:")
        for a in adapters:
            print(f"  {a}")
        return 0

    if cmd == "convert":
        if len(args) < 2:
            print("Usage: cfdx-io convert <case_path> [--output DIR] [--solver NAME]")
            return 1

        case_path = args[1]
        output_path = None
        solver = None

        i = 2
        while i < len(args):
            if args[i] == "--output" and i + 1 < len(args):
                output_path = args[i + 1]
                i += 2
            elif args[i] == "--solver" and i + 1 < len(args):
                solver = args[i + 1]
                i += 2
            else:
                i += 1

        detected = detect_adapter(case_path) if solver is None else solver
        print(f"Converting: {case_path}")
        print(f"Detected solver: {detected}")

        try:
            result = convert(case_path, output_path=output_path, solver=solver)
        except ValueError as e:
            print(f"Error: {e}")
            return 1

        if result.gap_report.has_blocking():
            print("\nConversion BLOCKED due to blocking incompatibilities.")
            print("See gap_analysis.md for details.")
            return 1

        print("\nConversion completed successfully.")
        print(f"  Supported features: {result.gap_report.n_supported()}")
        print(f"  Approximated: {result.gap_report.n_approximated()}")
        print(f"  Unsupported (non-blocking): {result.gap_report.n_unsupported_nonblocking()}")
        print(f"  Blocking issues: {result.gap_report.n_unsupported_blocking()}")

        if result.mesh:
            print(f"  Mesh: {result.mesh['points'].shape[0]} vertices")

        if output_path:
            print(f"\nReports written to: {output_path}")
            print(f"  - gap_analysis.md")
            print(f"  - gap_analysis.json")
            print(f"  - case_setup.json")
            print(f"  - mesh.vtk")

        return 0

    print(f"Unknown command: {cmd}")
    print("Commands: convert, list-adapters")
    return 1


if __name__ == "__main__":
    sys.exit(main())
