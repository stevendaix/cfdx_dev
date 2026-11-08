#!/usr/bin/env python3
"""Collect evidence without ever inferring physical qualification from CTest PASS."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/cfdx/python'))
from cfdx.qualification_evidence import qualification_evidence_manifest


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--registry', required=True)
    p.add_argument('--evidence-root', required=True)
    p.add_argument('--commit', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args(argv)
    output = Path(args.output)
    manifest = qualification_evidence_manifest(args.registry, args.evidence_root, args.commit)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
