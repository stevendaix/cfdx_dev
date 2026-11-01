#!/usr/bin/env python3
"""Validate independently produced CFDX NACA run summaries."""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

FIELDS = (
    "mesh",
    "strategy",
    "converged",
    "wall_seconds",
    "true_momentum_residual",
    "true_continuity_residual",
    "mass_imbalance",
    "cl",
    "cd",
    "cp_error",
    "cf_error",
)


def analyze(path):
    with open(path, newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not set(FIELDS).issubset(set(reader.fieldnames or [])):
            raise ValueError("Required columns: " + ", ".join(FIELDS))
        rows = list(reader)

    if not rows:
        raise ValueError("No runs found")

    output = []
    for row in rows:
        if row["mesh"] not in ("N64", "N128", "N256"):
            raise ValueError("Unknown mesh " + row["mesh"])
        if not row["strategy"].strip():
            raise ValueError("strategy must not be empty")
        if row["converged"].lower() not in ("true", "false"):
            raise ValueError("converged must be true or false")

        numeric = {}
        for key in FIELDS[3:]:
            value = float(row[key])
            if not math.isfinite(value):
                raise ValueError(f"{row['mesh']}: {key} nonfinite")
            numeric[key] = value

        nonnegative = (
            "wall_seconds",
            "true_momentum_residual",
            "true_continuity_residual",
            "mass_imbalance",
            "cp_error",
            "cf_error",
        )
        if any(numeric[key] < 0.0 for key in nonnegative):
            raise ValueError("Invalid negative measurement")

        output.append(
            {
                "mesh": row["mesh"],
                "strategy": row["strategy"],
                "converged": row["converged"].lower() == "true",
                **numeric,
            }
        )

    return {
        "results": output,
        "note": (
            "Input data reproduced only. Converged flags are user-supplied, "
            "not independently certified."
        ),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--out", default="naca_comparison.json")
    args = parser.parse_args()
    try:
        report = analyze(args.csv)
        Path(args.out).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"Validated {len(report['results'])} rows; wrote {args.out}")
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)
