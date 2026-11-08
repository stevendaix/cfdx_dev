"""Strict CI qualification evidence manifest; registry PASS is never promoted.

Each case has an independent proof record with the exact HEAD, tests and
SHA-256 of source artifacts. Missing or unmatched records remain UNVERIFIED.
"""
from __future__ import annotations

import hashlib
import math
import json
import re
from pathlib import Path

from .qualification_center import read_qualification_registry

_SHA = re.compile(r"[0-9a-f]{40}\Z")


def _hash(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def qualification_evidence_manifest(registry_path, evidence_root, commit_sha: str):
    """Only a verified per-case PASS certificate can change case proof state.

    A proof is a JSON file named by registry case id, with exact commit_sha,
    case_id, status PASS, and a nonempty list of test/evidence artifact names.
    Missing, stale, malformed and inconclusive proofs are never PASS.
    """
    if not _SHA.fullmatch(commit_sha):
        raise ValueError("required exact 40-character lowercase commit SHA")
    registry = read_qualification_registry(registry_path)
    root = Path(evidence_root).resolve()
    proofs = []
    for case in registry.cases:
        path = root / f"{case.id}.json"
        status, reason, evidence_hash = "UNVERIFIED", "no current-commit proof", None
        if path.is_file() and not path.is_symlink() and path.stat().st_size <= 256 * 1024:
            try:
                doc = json.loads(path.read_text(encoding="utf-8"))
                if (doc.get("case_id") == case.id and doc.get("commit_sha") == commit_sha and
                    doc.get("status") == "PASS" and doc.get("qualifier") == "physical" and
                    isinstance(doc.get("artifacts"), list) and doc["artifacts"]):
                    names = doc["artifacts"]
                    if any(not isinstance(v, str) or not v or Path(v).name != v for v in names):
                        raise ValueError("unsafe evidence path")
                    hashes = doc.get("artifact_sha256")
                    checks = doc.get("checks")
                    if not isinstance(hashes, dict) or set(hashes) != set(names):
                        raise ValueError("missing artifact hashes")
                    if not isinstance(checks, list) or not checks:
                        raise ValueError("missing physical oracle checks")
                    for check in checks:
                        if not isinstance(check, dict) or not isinstance(check.get("name"), str):
                            raise ValueError("malformed physical check")
                        measured, limit = check.get("error"), check.get("limit")
                        if any(isinstance(x, bool) or not isinstance(x, (float, int)) or
                               not math.isfinite(x) or x < 0 for x in (measured, limit)):
                            raise ValueError("physical check requires finite nonnegative error and limit")
                        if measured > limit or check.get("passed") is not True:
                            raise ValueError("physical oracle failed")
                    for name in names:
                        artifact = root / name
                        if not artifact.is_file() or artifact.is_symlink() or artifact.stat().st_size > 512 * 1024 * 1024:
                            raise ValueError("missing or oversized evidence artifact")
                        if _hash(artifact) != hashes[name]:
                            raise ValueError("physical evidence checksum mismatch")
                    status, reason, evidence_hash = "PASS_EVIDENCED", "current HEAD + checked physical artifacts", _hash(path)
                else:
                    reason = "proof incomplete, not physical, or on a different HEAD"
            except (OSError, ValueError, TypeError, AttributeError):
                reason = "malformed or incomplete proof"
        proofs.append({"case_id": case.id, "declared_status": case.status,
                       "proof_status": status, "reason": reason, "evidence_sha256": evidence_hash})
    return {"format": "CFDX-CI-QUALIFICATION-EVIDENCE", "schema_version": 1,
            "commit_sha": commit_sha, "cases": proofs,
            "physically_evidenced": sum(p["proof_status"] == "PASS_EVIDENCED" for p in proofs),
            "total": len(proofs)}
