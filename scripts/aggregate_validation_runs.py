"""Aggregate validation CSV files deterministically with provenance checks."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if (
            not reader.fieldnames
            or any(not name or not name.strip() for name in reader.fieldnames)
        ):
            raise ValueError(f"{path}: missing or invalid CSV header")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"{path}: duplicate CSV header")
        if "source_file" in reader.fieldnames:
            raise ValueError(
                f"{path}: input cannot override reserved source_file provenance"
            )
        rows = list(reader)
        if any(
            None in row or any(value is None for value in row.values()) for row in rows
        ):
            raise ValueError(f"{path}: row does not match the CSV header")
        return list(reader.fieldnames), rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--pattern", default="*.csv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--key", nargs="+", default=["case"])
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {args.input_dir}")

    output = args.output.resolve()
    files = sorted(
        p for p in args.input_dir.glob(args.pattern)
        if p.is_file() and p.resolve() != output
    )
    if not files:
        raise SystemExit(f"No files matching {args.pattern} in {args.input_dir}")

    rows: list[dict[str, str]] = []
    fields: list[str] = []
    seen: set[tuple[str, ...]] = set()

    for path in files:
        file_fields, file_rows = read_rows(path)
        for field in file_fields:
            if field not in fields:
                fields.append(field)
        for row in file_rows:
            missing = [key for key in args.key if key not in row]
            if missing:
                raise ValueError(
                    f"{path}: missing key columns: {', '.join(missing)}"
                )
            key = tuple(row[key].strip() for key in args.key)
            if any(not part for part in key):
                raise ValueError(f"{path}: validation key contains an empty value")
            if key in seen:
                raise ValueError(
                    f"duplicate validation key {key} found while reading {path}"
                )
            seen.add(key)
            row = dict(row)
            row["source_file"] = str(path.relative_to(args.input_dir))
            if "source_file" not in fields:
                fields.append("source_file")
            rows.append(row)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
