"""Read-only qualification registry contract for the CFDX Workbench.

The registry is a *declaration* of engineering progress, not proof that a
particular solver checkout has passed a validation campaign.  This module
never promotes a case to PASS or changes a scientific acceptance threshold.
It is Qt-free and can be used by CI, CLI and the desktop application.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

REGISTRY_FILENAME = "CFDX_QUALIFICATION_REGISTRY.json"
REGISTRY_RELATIVE_PATH = Path("docs/validation") / REGISTRY_FILENAME
MAX_REGISTRY_BYTES = 4 * 1024 * 1024
MAX_REGISTRY_CASES = 5_000
MAX_FIELD_CHARACTERS = 8_192
SUPPORTED_SCHEMA_VERSION = "1.0"
CASE_ID_PATTERN = re.compile(r"[A-Z0-9][A-Z0-9_-]{1,79}\Z")
STATUSES = ("PLANNED", "READY", "RUNNING", "DIAGNOSTIC", "PASS", "BLOCKED")


@dataclass(frozen=True, slots=True)
class QualificationCase:
    id: str
    domain: str
    name: str
    status: str
    details: tuple[tuple[str, str], ...]

    def matches(self, query: str = "", *, domain: str = "", status: str = "") -> bool:
        if domain and self.domain != domain:
            return False
        if status and self.status != status:
            return False
        needle = query.casefold().strip()
        return not needle or any(
            needle in value.casefold()
            for value in (self.id, self.domain, self.name, self.status)
        )


@dataclass(frozen=True, slots=True)
class QualificationRegistry:
    source: Path
    campaign: str
    cases: tuple[QualificationCase, ...]
    policy: tuple[tuple[str, str], ...]
    source_signature: tuple[int, int, int, int]

    def filter(
        self, query: str = "", *, domain: str = "", status: str = ""
    ) -> tuple[QualificationCase, ...]:
        return tuple(
            case
            for case in self.cases
            if case.matches(query, domain=domain, status=status)
        )

    def counts(self) -> dict[str, int]:
        return {
            status: sum(case.status == status for case in self.cases)
            for status in STATUSES
        }

    def policy_for(self, status: str) -> str:
        return dict(self.policy)[status]


def _unique_json_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject ambiguous duplicate JSON keys rather than trusting the last one."""
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str):
    raise ValueError(f"non-finite JSON constant: {value}")


def _bounded_string(value: object, label: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str) or len(value) > MAX_FIELD_CHARACTERS:
        raise ValueError(f"{label} must be a bounded string")
    if not allow_empty and not value.strip():
        raise ValueError(f"{label} must not be empty")
    return value


def _signature(source: Path) -> tuple[int, int, int, int]:
    stat = source.stat()
    return int(stat.st_dev), int(stat.st_ino), int(stat.st_size), int(stat.st_mtime_ns)


def read_qualification_registry(source: str | Path) -> QualificationRegistry:
    """Parse a bounded registry snapshot without altering its source.

    A changed file is rejected instead of returning a mixed/stale snapshot.
    A declared PASS is *not* evidence of a PASS on the running Git HEAD.
    """
    path = Path(source).resolve()
    before = _signature(path)
    if before[2] > MAX_REGISTRY_BYTES:
        raise ValueError("qualification registry exceeds the 4 MiB input limit")
    payload = path.read_bytes()
    if len(payload) != before[2] or _signature(path) != before:
        raise RuntimeError("qualification registry changed during reading")
    raw = json.loads(
        payload.decode("utf-8"),
        object_pairs_hook=_unique_json_keys,
        parse_constant=_reject_json_constant,
    )
    if (
        not isinstance(raw, dict)
        or raw.get("schema_version") != SUPPORTED_SCHEMA_VERSION
    ):
        raise ValueError("unsupported qualification registry schema version")
    campaign = _bounded_string(raw.get("campaign"), "campaign")
    policies = raw.get("status_policy")
    if not isinstance(policies, dict) or set(policies) != set(STATUSES):
        raise ValueError("qualification registry has incomplete status policy")
    policy = tuple(
        (status, _bounded_string(policies[status], f"status_policy.{status}"))
        for status in STATUSES
    )
    data = raw.get("cases")
    if not isinstance(data, list) or len(data) > MAX_REGISTRY_CASES:
        raise ValueError("invalid qualification registry case list")
    cases: list[QualificationCase] = []
    seen: set[str] = set()
    for index, item in enumerate(data):
        if not isinstance(item, dict):
            raise ValueError(f"case[{index}] must be an object")  # noqa: TRY004
        case_id = _bounded_string(item.get("id"), f"case[{index}].id")
        if not CASE_ID_PATTERN.fullmatch(case_id):
            raise ValueError(f"invalid qualification case id: {case_id!r}")
        if case_id in seen:
            raise ValueError(f"duplicate qualification case id: {case_id!r}")
        seen.add(case_id)
        domain = _bounded_string(item.get("domain"), f"case[{index}].domain")
        name = _bounded_string(item.get("case"), f"case[{index}].case")
        status = _bounded_string(item.get("status"), f"case[{index}].status")
        if status not in STATUSES:
            raise ValueError(
                f"unsupported qualification status for {case_id}: {status}"
            )
        details = tuple(
            (key, _bounded_string(value, f"{case_id}.{key}", allow_empty=True))
            for key, value in item.items()
            if key not in {"id", "domain", "case", "status"}
        )
        cases.append(QualificationCase(case_id, domain, name, status, details))
    return QualificationRegistry(path, campaign, tuple(cases), policy, before)


def discover_qualification_registry(start: str | Path | None = None) -> Path | None:
    """Find the existing repository registry; do not generate a second copy."""
    roots = (
        [Path(start).resolve()]
        if start is not None
        else [Path.cwd().resolve(), Path(__file__).resolve()]
    )
    for root in roots:
        for parent in (root, *root.parents):
            candidate = parent / REGISTRY_RELATIVE_PATH
            if candidate.is_file():
                return candidate
    return None
