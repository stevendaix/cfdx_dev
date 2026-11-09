"""Read-only cross-case qualification snapshot, with explicit provenance limits.

Example:
  PYTHONPATH=src/cfdx/python python scripts/cfdx_qualification_status.py --status DIAGNOSTIC
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# This script is also usable from an uninstalled source checkout.
_PACKAGE_ROOT = Path(__file__).resolve().parents[1] / "src/cfdx/python"
# The repository also has scripts/cfdx.py; it would shadow the package when
# the source directory is present on PYTHONPATH but ordered after scripts/.
if str(_PACKAGE_ROOT) in sys.path:
    sys.path.remove(str(_PACKAGE_ROOT))
sys.path.insert(0, str(_PACKAGE_ROOT))

from cfdx.qualification_center import (
    STATUSES,
    discover_qualification_registry,
    read_qualification_registry,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, help="registry JSON to inspect")
    parser.add_argument("--status", choices=STATUSES, default="")
    parser.add_argument("--domain", default="", help="exact domain filter")
    parser.add_argument("--query", default="", help="case-insensitive ID/name search")
    parser.add_argument("--format", choices=("table", "json"), default="table")
    args = parser.parse_args(argv)
    source = args.registry or discover_qualification_registry()
    if source is None:
        parser.error("qualification registry not found; pass --registry")
    try:
        registry = read_qualification_registry(source)
    except (OSError, RuntimeError, ValueError, UnicodeError) as exc:
        parser.exit(2, f"qualification registry rejected: {exc}\n")
    cases = registry.filter(args.query, domain=args.domain, status=args.status)
    if args.format == "json":
        output = {
            "source": str(registry.source),
            "campaign": registry.campaign,
            "status_semantics": "declarations_only_not_certified_on_current_head",
            "total_cases": len(registry.cases),
            "counts": registry.counts(),
            "filtered_cases": [
                {
                    "id": case.id,
                    "case": case.name,
                    "domain": case.domain,
                    "declared_status": case.status,
                    "head_evidence": "unverified",
                }
                for case in cases
            ],
        }
        print(json.dumps(output, ensure_ascii=False, indent=2))
    else:
        print(f"CFDX Qualification Center — {len(cases)}/{len(registry.cases)} cases")
        print("Registry declarations ONLY; current-HEAD evidence is NOT verified.")
        print(f"Registry: {registry.source}")
        print("  ".join(f"{name}={count}" for name, count in registry.counts().items()))
        for case in cases:
            print(f"  {case.id:<32} {case.status:<12} {case.domain} — {case.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
