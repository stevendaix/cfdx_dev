#!/usr/bin/env python3
"""CFDX Convert CLI - Python-side entry point for format conversion.

This module provides the Python-side conversion logic called by the C++ CLI.
It uses the unified converter API from cfdx.io.pipeline and writes the output
as a .cfdx.h5 file with proper schema.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from cfdx.io.pipeline import convert, detect_adapter, get_available_adapters
from cfdx.io.interfaces import ConversionResult
from cfdx.io.schema import CFDX_SCHEMA_VERSION, Severity


def write_cfdx_h5(result: ConversionResult, output_path: Path) -> None:
    """Write ConversionResult to a CFDX HDF5 file."""
    with h5py.File(output_path, "w") as f:
        f.attrs["cfdx_version"] = CFDX_SCHEMA_VERSION
        f.attrs["source_format"] = result.source.solver if result.source.solver else "unknown"
        f.attrs["source_version"] = result.source.version if result.source.version else ""
        f.attrs["source_file"] = result.source.case_path if result.source.case_path else ""

        if result.mesh is not None and result.mesh.get("points") is not None:
            points = result.mesh["points"]
            cells = result.mesh.get("cells", {})

            mesh_grp = f.create_group("mesh/topology")
            mesh_grp.attrs["cell_count"] = len(cells) if cells else 0
            mesh_grp.attrs["point_count"] = points.shape[0] if points is not None else 0
            mesh_grp.attrs["dimension"] = points.shape[1] if points is not None and len(points.shape) > 1 else 3

            pts_ds = mesh_grp.create_dataset("points", data=points, compression="gzip")
            pts_ds.attrs["description"] = "Mesh vertex coordinates (x, y, z)"

            if cells:
                conn_grp = mesh_grp.create_group("connectivity")
                for cell_type, conn in cells.items():
                    conn_grp.create_dataset(cell_type, data=conn, compression="gzip")

        if result.scalar_fields:
            fields_grp = f.create_group("fields/scalar")
            for name, data in result.scalar_fields:
                fields_grp.create_dataset(name, data=data, compression="gzip")

        if result.vec_fields:
            vec_grp = f.create_group("fields/vector")
            for name, data in result.vec_fields:
                vec_grp.create_dataset(name, data=data, compression="gzip")

        if result.setup is not None:
            setup_grp = f.create_group("case_setup")
            setup_dict = result.setup.model_dump(mode="json")
            setup_grp.attrs["json"] = json.dumps(setup_dict)

        if result.gap_report is not None:
            gap_grp = f.create_group("gap_analysis")
            gap_grp.attrs["json"] = json.dumps({
                "supported": result.gap_report.n_supported(),
                "approximated": result.gap_report.n_approximated(),
                "unsupported_nonblocking": result.gap_report.n_unsupported_nonblocking(),
                "unsupported_blocking": result.gap_report.n_unsupported_blocking(),
                "findings": [
                    {
                        "feature": f.feature,
                        "severity": f.severity.value,
                        "message": f.detail,
                        "suggestion": f.suggestion
                    }
                    for f in result.gap_report.findings
                ]
            })

        if result.available_fields:
            avail_grp = f.create_group("available_fields")
            for category, fields in result.available_fields.items():
                avail_grp.attrs[category] = json.dumps(fields)

        if result.field_data_source:
            f.attrs["field_data_source"] = result.field_data_source


def build_summary(result: ConversionResult, output_path: Path) -> dict:
    """Build a JSON-serializable summary of the conversion."""
    summary = {
        "source_format": result.source.solver if result.source.solver else "unknown",
        "source_version": result.source.version if result.source.version else "unknown",
        "source_file": str(result.source.case_path) if result.source.case_path else "unknown",
        "output_file": str(output_path),
        "success": result.success,
        "mesh": {},
        "fields": {
            "scalar": len(result.scalar_fields),
            "vector": len(result.vec_fields)
        },
        "gap_analysis": {
            "supported": result.gap_report.n_supported() if result.gap_report else 0,
            "approximated": result.gap_report.n_approximated() if result.gap_report else 0,
            "unsupported_nonblocking": result.gap_report.n_unsupported_nonblocking() if result.gap_report else 0,
            "unsupported_blocking": result.gap_report.n_unsupported_blocking() if result.gap_report else 0,
        }
    }

    if result.mesh is not None and result.mesh.get("points") is not None:
        points = result.mesh["points"]
        cells = result.mesh.get("cells", {})
        total_cells = sum(len(v) for v in cells.values()) if cells else 0
        summary["mesh"] = {
            "points": points.shape[0] if points is not None else 0,
            "cells": total_cells,
            "dimension": points.shape[1] if points is not None and len(points.shape) > 1 else 3,
            "cell_types": list(cells.keys()) if cells else []
        }

    if result.scalar_fields:
        summary["fields"]["scalar_names"] = [name for name, _ in result.scalar_fields]
    if result.vec_fields:
        summary["fields"]["vector_names"] = [name for name, _ in result.vec_fields]

    if result.setup is not None:
        summary["case_setup"] = {
            "boundary_patches": len(result.setup.boundary_conditions) if result.setup.boundary_conditions else 0,
            "materials": len(result.setup.materials) if result.setup.materials else 0,
            "initial_conditions": 1 if result.setup.initial_condition else 0,
        }

    return summary


def convert_cli(source: str, output: str) -> int:
    """Main conversion function called by C++ bridge."""
    source_path = Path(source)
    output_path = Path(output)

    print(f"Converting {source_path} -> {output_path}", file=sys.stderr)

    try:
        result = convert(source_path, output_path=output_path.parent if output_path.parent != Path(".") else None)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unexpected error during conversion: {e}", file=sys.stderr)
        return 1

    if result.gap_report.has_blocking():
        print("Conversion has blocking incompatibilities:", file=sys.stderr)
        for finding in result.gap_report.findings:
            if finding.severity == Severity.UNSUPPORTED_BLOCK:
                print(f"  BLOCKING: {finding.feature} - {finding.detail}", file=sys.stderr)
        return 1

    try:
        write_cfdx_h5(result, output_path)
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        return 1

    summary = build_summary(result, output_path)
    print(json.dumps(summary, indent=2))

    return 0


def normalize_cli(source: str, output: str) -> int:
    """Normalize/re-export a CFDX file."""
    source_path = Path(source)
    output_path = Path(output)

    print(f"Normalizing {source_path} -> {output_path}", file=sys.stderr)

    try:
        with h5py.File(source_path, "r") as src:
            with h5py.File(output_path, "w") as dst:
                for attr_name, attr_val in src.attrs.items():
                    dst.attrs[attr_name] = attr_val

                def copy_group(name, obj):
                    if isinstance(obj, h5py.Group):
                        dst.create_group(name)
                    elif isinstance(obj, h5py.Dataset):
                        src.copy(name, dst, name=name)

                src.visititems(copy_group)

        summary = {
            "source_format": "CFDX",
            "source_file": str(source_path),
            "output_file": str(output_path),
            "success": True,
            "action": "normalized"
        }
        print(json.dumps(summary, indent=2))
        return 0
    except Exception as e:
        print(f"Error normalizing CFDX file: {e}", file=sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]

    if not args:
        print("Usage: cfdx-convert <command> [args...]", file=sys.stderr)
        print("Commands:", file=sys.stderr)
        print("  convert <source> <output>  - Convert from external format to CFDX", file=sys.stderr)
        print("  normalize <source> <output> - Normalize/re-export CFDX file", file=sys.stderr)
        print("  list-adapters              - List available solver adapters", file=sys.stderr)
        return 1

    cmd = args[0]

    if cmd == "list-adapters":
        adapters = get_available_adapters()
        print("Available solver adapters:")
        for a in adapters:
            print(f"  {a}")
        return 0

    if cmd == "convert":
        if len(args) != 3:
            print("Usage: cfdx-convert convert <source> <output>", file=sys.stderr)
            return 1

        source = args[1]
        output = args[2]

        if source.endswith(".cfdx.h5") or Path(source).name == ".cfdx.h5":
            return normalize_cli(source, output)

        return convert_cli(source, output)

    if cmd == "normalize":
        if len(args) != 3:
            print("Usage: cfdx-convert normalize <source> <output>", file=sys.stderr)
            return 1
        return normalize_cli(args[1], args[2])

    print(f"Unknown command: {cmd}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())