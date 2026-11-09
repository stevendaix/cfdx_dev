"""All-case qualification portfolio: metadata, CTest runs, and physical evidence.

This is deliberately a *tracker* rather than a physical-certification engine.
A CTest pass cannot promote a registry case. The authoritative independent
physical-evidence verifier is scripts/qualification_evidence_gate.py.

The tracker is read-only with respect to the registry, physics and tolerances.
Absent artifacts, skipped tests, duplicate JUnit records and stale reports
never turn into PASS. A source archive hash describes source provenance, not
proof that the executed binary was rebuilt from that source.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import qualification_evidence_gate as evidence_gate

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "docs/validation/CFDX_QUALIFICATION_REGISTRY.json"
LINKS = ROOT / "docs/validation/QUALIFICATION_PORTFOLIO_LINKS.json"
SCHEMA = "cfdx-qualification-portfolio/v1"


class PortfolioError(ValueError):
    """Input sources or the explicit qualification mapping are inconsistent."""


def read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise PortfolioError(f"expected JSON object: {path}")
    return obj


def validate_links(registry: dict, mapping: dict) -> None:
    cases = registry.get("cases")
    links = mapping.get("cases")
    if mapping.get("schema_version") != 1 or not isinstance(links, dict):
        raise PortfolioError("invalid case-link catalog schema")
    if not isinstance(cases, list) or not cases:
        raise PortfolioError("empty registry")
    ids = [entry.get("id") for entry in cases]
    if any(not isinstance(i, str) for i in ids) or len(ids) != len(set(ids)):
        raise PortfolioError("registry case IDs missing or duplicated")
    if set(ids) != set(evidence_gate.REQUIRED_QOI):
        raise PortfolioError(
            "physical-evidence verifier no longer covers every case ID"
        )
    if set(ids) != set(links):
        raise PortfolioError(
            f"case mapping drift: missing={sorted(set(ids) - set(links))}; "
            f"extra={sorted(set(links) - set(ids))}"
        )
    for case_id, link in links.items():
        if not isinstance(link, dict):
            raise PortfolioError(f"invalid mapping for {case_id}")
        for field in ("contract_tests", "campaign_tests", "reference_workflows"):
            items = link.get(field)
            if (
                not isinstance(items, list)
                or any(not isinstance(n, str) or not n for n in items)
                or len(items) != len(set(items))
            ):
                raise PortfolioError(f"invalid {field} for {case_id}")
        if set(link["contract_tests"]) & set(link["campaign_tests"]):
            raise PortfolioError(f"duplicate contract/campaign test for {case_id}")
        if (
            not link["contract_tests"]
            or not isinstance(link.get("next_action"), str)
            or not link["next_action"].strip()
        ):
            raise PortfolioError(f"missing contract or action for {case_id}")
        if not isinstance(link.get("domain_group"), str) or not link["domain_group"]:
            raise PortfolioError(f"missing domain_group for {case_id}")
        for workflow in link["reference_workflows"]:
            if (
                Path(workflow).name != workflow
                or not workflow.endswith(".yml")
                or not (ROOT / ".github/workflows" / workflow).is_file()
            ):
                raise PortfolioError(
                    f"invalid or missing reference workflow {workflow} for {case_id}"
                )
    backlog = mapping.get("unregistered_backlog")
    if not isinstance(backlog, list):
        raise PortfolioError("unregistered_backlog must be a list")
    extras = [row.get("id") for row in backlog if isinstance(row, dict)]
    if (
        len(extras) != len(backlog)
        or len(extras) != len(set(extras))
        or set(extras) & set(ids)
    ):
        raise PortfolioError("backlog IDs missing, duplicated, or registered already")
    for item in backlog:
        if not isinstance(item.get("reason"), str) or not item["reason"]:
            raise PortfolioError(f"backlog missing reason: {item!r}")
        if not isinstance(item.get("existing_checks"), list):
            raise PortfolioError(f"backlog missing checks: {item!r}")


def read_manifest(path: Path | None) -> dict[str, bool] | None:
    if path is None:
        return None
    doc = read_json(path)
    tests = doc.get("tests")
    if not isinstance(tests, list) or not tests:
        raise PortfolioError("CTest JSON manifest is empty or invalid")
    names = [item.get("name") for item in tests if isinstance(item, dict)]
    if (
        len(names) != len(tests)
        or any(not isinstance(n, str) or not n for n in names)
        or len(names) != len(set(names))
    ):
        raise PortfolioError("CTest manifest has empty or duplicated test names")
    # CTest accepts DISABLED tests in --show-only. They must not be counted as
    # executable or even a valid mandatory anchor in the strict inventory.
    enabled: dict[str, bool] = {}
    for item in tests:
        properties = {p.get("name"): p.get("value") for p in item.get("properties", [])}
        enabled[item["name"]] = properties.get("DISABLED") not in (
            True,
            "TRUE",
            "ON",
            "1",
            1,
        )
    return enabled


def read_junits(paths: list[Path]) -> dict[str, str]:
    results: dict[str, str] = {}
    for path in paths:
        root = ET.parse(path).getroot()
        if root.tag not in ("testsuite", "testsuites"):
            raise PortfolioError(f"invalid JUnit XML root in {path}")
        for node in root.iter("testcase"):
            name = node.get("name")
            if not name:
                raise PortfolioError(f"unnamed JUnit testcase in {path}")
            if node.find("failure") is not None or node.find("error") is not None:
                state = "FAIL"
            elif node.find("skipped") is not None or node.get("status") in (
                "notrun",
                "skipped",
            ):
                state = "SKIPPED"
            else:
                state = "PASS"
            if name in results:
                results[name] = "DUPLICATE"
            else:
                results[name] = state
    return results


def test_state(
    name: str, registered: dict[str, bool] | set[str] | None, junit: dict[str, str]
) -> str:
    if registered is not None and name not in registered:
        return "NOT_REGISTERED"
    if isinstance(registered, dict) and registered.get(name) is False:
        return "DISABLED"
    return junit.get(name, "NOT_RUN")


def aggregate_status(records: dict[str, str]) -> str:
    if not records:
        return "NO_DEDICATED_TEST"
    all_statuses = set(records.values())
    if "DUPLICATE" in all_statuses:
        return "AMBIGUOUS"
    if "FAIL" in all_statuses:
        return "FAIL"
    if "NOT_REGISTERED" in all_statuses:
        return "NOT_REGISTERED"
    if "DISABLED" in all_statuses:
        return "DISABLED"
    if "SKIPPED" in all_statuses:
        return "SKIPPED"
    if all_statuses == {"PASS"}:
        return "PASS"
    if all_statuses == {"NOT_RUN"}:
        return "NOT_RUN"
    return "PARTIAL"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def make_portfolio(
    registry: dict,
    mapping: dict,
    registered: dict[str, bool] | set[str] | None,
    junit: dict[str, str],
    evidence_dir: Path,
    source_archive: Path | None = None,
    source_revision: str | None = None,
) -> tuple[dict, dict]:
    validate_links(registry, mapping)
    evidence = evidence_gate.evaluate_registry(registry, evidence_dir)
    ev_by_id = {row["id"]: row for row in evidence["cases"]}
    rows: list[dict[str, Any]] = []
    for case in registry["cases"]:
        ident = case["id"]
        linked = mapping["cases"][ident]
        contract = {
            n: test_state(n, registered, junit) for n in linked["contract_tests"]
        }
        campaign = {
            n: test_state(n, registered, junit) for n in linked["campaign_tests"]
        }
        physical = ev_by_id[ident]
        physical_ok = physical["verification_status"] == "EVIDENCE_GATES_PASS"
        missing_checks = [
            item["check"] for item in physical["checks"] if item["status"] != "PASS"
        ]
        rows.append(
            {
                "id": ident,
                "domain_group": linked["domain_group"],
                "name": case.get("case", ident),
                "registry_status": case.get("status", "UNKNOWN"),
                "contracts": contract,
                "contract_status": aggregate_status(contract),
                "campaigns": campaign,
                "campaign_status": aggregate_status(campaign),
                "physical_evidence_status": physical["verification_status"],
                "physical_failed_checks": missing_checks,
                "reference_workflows": linked["reference_workflows"],
                "next_action": linked["next_action"],
                "tracking_note": (
                    "Physical evidence gates passed; independent physical review and exact-source correspondence still required"
                    if physical_ok
                    else "Not physically qualified on supplied evidence"
                ),
            }
        )
    backlog = []
    for item in mapping["unregistered_backlog"]:
        checks = {
            name: test_state(name, registered, junit)
            for name in item["existing_checks"]
        }
        backlog.append(
            {
                "id": item["id"],
                "domain_group": item["domain"],
                "registry_status": "NOT_IN_REGISTRY",
                "checks": checks,
                "checks_status": aggregate_status(checks),
                "reason": item["reason"],
                "physical_evidence_status": "NOT_APPLICABLE_UNREGISTERED",
            }
        )
    missed = sorted(
        {
            name
            for row in rows
            for name, result in {**row["contracts"], **row["campaigns"]}.items()
            if result in ("NOT_REGISTERED", "DISABLED")
        }
    )
    by_domain = {}
    for domain in sorted({row["domain_group"] for row in rows}):
        subset = [row for row in rows if row["domain_group"] == domain]
        by_domain[domain] = {
            "cases": len(subset),
            "registry_legacy_pass": sum(r["registry_status"] == "PASS" for r in subset),
            "contract_run_pass": sum(r["contract_status"] == "PASS" for r in subset),
            "campaign_run_pass": sum(r["campaign_status"] == "PASS" for r in subset),
            "evidence_gates_pass": sum(
                r["physical_evidence_status"] == "EVIDENCE_GATES_PASS" for r in subset
            ),
        }
    board = {
        "schema": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "provenance": {
            "source_archive": (
                str(source_archive) if source_archive is not None else None
            ),
            "source_archive_sha256": (
                file_sha256(source_archive) if source_archive is not None else None
            ),
            "declared_source_revision": source_revision,
            "physical_evidence_directory": str(evidence_dir),
            "ctest_manifest_provided": registered is not None,
            "junit_files_provided": len(junit) > 0,
            "binary_correspondence_proven": False,
        },
        "summary": {
            "registered_cases": len(rows),
            "registry_status": dict(
                sorted(Counter(row["registry_status"] for row in rows).items())
            ),
            "physical_evidence_gates_pass": evidence["coverage"]["evidence_pass"],
            "registered_missing_ctest_anchors": missed,
            "unregistered_backlog": len(backlog),
            "ctest_inventory_size": len(registered) if registered is not None else None,
            "provided_junit_test_results": len(junit),
            "ctest_junit_pass": sum(s == "PASS" for s in junit.values()),
            "ctest_junit_fail_or_ambiguous": sum(
                s in ("FAIL", "SKIPPED", "DUPLICATE") for s in junit.values()
            ),
            "by_domain": by_domain,
        },
        "cases": rows,
        "unregistered_backlog": backlog,
        "warnings": [
            "Legacy registry PASS is historical: not a qualification on the provided source.",
            "CTest PASS proves only the executed test's claims; it does not certify a CFD model.",
            "This tracker never promotes physical status. Check evidence gate, exact source/build correspondence and independent validation manually.",
            "Reference-only workflows do not prove that CFDX solver cases have executed.",
        ],
    }
    return board, evidence


def markdown_table_cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def markdown_report(report: dict, previous: dict | None = None) -> str:
    summary = report["summary"]
    src = report["provenance"]
    out = [
        "# CFDX — Qualification simultanée : tableau de bord",
        "",
        f"Généré : {report['generated_at_utc']}.",
        (
            f"Source archive SHA-256 : `{src['source_archive_sha256'] or 'NON FOURNIE'}`. "
            "Correspondance binaire/source : **non attestée**."
        ),
        f"Révision source déclarée : `{markdown_table_cell(src['declared_source_revision'] or 'NON FOURNIE')}`.",
        "",
        (
            f"**{summary['registered_cases']} cas enregistrés ; {summary['unregistered_backlog']} "
            "campagnes supplémentaires à inscrire.** "
            f"Preuves physiques vérifiables dans le dossier fourni : **{summary['physical_evidence_gates_pass']}/{summary['registered_cases']}**. "
            f"Manifest CTest : **{summary['ctest_inventory_size'] if summary['ctest_inventory_size'] is not None else 'non fourni'}** tests, "
            f"JUnit importés : **{summary['provided_junit_test_results']}**."
        ),
        "",
        (
            "> `PASS` de la colonne registre = historique. `PASS` de test = résultat JUnit. "
            "Les deux sont indépendants de `EVIDENCE_GATES_PASS`, et ce dernier n'est pas une certification scientifique automatique."
        ),
        "",
        "## Vue par domaine",
        "",
        "| Domaine | Cas | PASS registre historique | Contrats tous exécutés PASS | Campagnes tous exécutées PASS | Preuves physiques |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for dom, row in summary["by_domain"].items():
        out.append(
            f"| {markdown_table_cell(dom)} | {row['cases']} | {row['registry_legacy_pass']} | "
            f"{row['contract_run_pass']} | {row['campaign_run_pass']} | {row['evidence_gates_pass']} |"
        )
    out += [
        "",
        "## Suivi détaillé des 25 cas",
        "",
        "| Cas | Domaine | Registre | Contrats CTest | Campagne CFDX CTest | Preuve physique | Prochaine action |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in report["cases"]:
        out.append(
            "| "
            + " | ".join(
                markdown_table_cell(v)
                for v in (
                    f"`{row['id']}`",
                    row["domain_group"],
                    row["registry_status"],
                    row["contract_status"],
                    row["campaign_status"],
                    row["physical_evidence_status"],
                    row["next_action"],
                )
            )
            + " |"
        )
    out += ["", "## Tests réellement fournis / non exécutés", ""]
    for row in report["cases"]:
        line = [f"- **`{row['id']}`** :"]
        line.append(
            "contrats "
            + (", ".join(f"`{n}`={s}" for n, s in row["contracts"].items()) or "aucun")
        )
        line.append(
            "campagnes "
            + (", ".join(f"`{n}`={s}" for n, s in row["campaigns"].items()) or "aucune")
        )
        if row["reference_workflows"]:
            line.append(
                "workflows référence uniquement "
                + ", ".join(f"`{w}`" for w in row["reference_workflows"])
            )
        out.append(" ".join(line))
    out += [
        "",
        "## Campagnes à créer / déclarer dans le registre",
        "",
        "| Proposition | Tests déjà présents | Manque |",
        "|---|---|---|",
    ]
    for row in report["unregistered_backlog"]:
        out.append(
            f"| `{row['id']}` | {row['checks_status']} | {markdown_table_cell(row['reason'])} |"
        )
    missed = summary["registered_missing_ctest_anchors"]
    out += [
        "",
        "## Contrat de suivi et de promotion",
        "",
        f"- Ancrages CTest déclarés, absents ou désactivés : **{len(missed)}**"
        + (f" ({', '.join(missed)})" if missed else "")
        + ".",
        "- Le dossier `qualification_evidence/<CAS>/` doit contenir des preuves réelles d'au moins trois maillages, la configuration, le binaire, le journal, la référence et les oracles hashés.",
        "- Chaque campagne en échec ou non exécutée reste visible : jamais de `PASS` par absence de résultat ni par héritage de la CI.",
        "- La comparaison aux références physiques et les bilans massiques/momentum/énergie/rayonnement sont obligatoires selon le cas.",
        "- Relancer le même tableau après chaque groupe de campagnes sur le même état du code; les tests de référence seuls ne changent pas les preuves du solveur.",
    ]
    if previous is not None:
        prev = {r["id"]: r for r in previous["cases"]}
        out += ["", "## Changements depuis le précédent bilan", ""]
        changed = []
        for cur in report["cases"]:
            before = prev.get(cur["id"])
            if before is None:
                changed.append(f"- `{cur['id']}` : nouveau cas")
                continue
            for key in (
                "registry_status",
                "contract_status",
                "campaign_status",
                "physical_evidence_status",
            ):
                if before.get(key) != cur.get(key):
                    changed.append(
                        f"- `{cur['id']}` {key}: `{before.get(key)}` → `{cur.get(key)}`"
                    )
        out += (
            changed
            if changed
            else ["Aucun changement de statut mesuré entre ces deux bilans."]
        )
    out += ["", "*Document de suivi, non certificat de validation physique.*", ""]
    return "\n".join(out)


def atomic_write(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(payload, encoding="utf-8")
    os.replace(tmp, path)


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "case_id",
                "domain",
                "registry_status",
                "contract_status",
                "campaign_status",
                "physical_evidence_status",
                "next_action",
            )
        )
        for row in rows:
            writer.writerow(
                (
                    row["id"],
                    row["domain_group"],
                    row["registry_status"],
                    row["contract_status"],
                    row["campaign_status"],
                    row["physical_evidence_status"],
                    row["next_action"],
                )
            )
    os.replace(tmp, path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=REGISTRY)
    parser.add_argument("--links", type=Path, default=LINKS)
    parser.add_argument(
        "--ctest-manifest", type=Path, help="ctest --show-only=json-v1 output"
    )
    parser.add_argument(
        "--junit",
        type=Path,
        action="append",
        default=[],
        help="actual CTest JUnit XML; can repeat",
    )
    parser.add_argument(
        "--evidence-dir", type=Path, default=ROOT / "qualification_evidence"
    )
    parser.add_argument(
        "--source-archive", type=Path, help="optional immutable source ZIP, sha256 only"
    )
    parser.add_argument(
        "--source-revision",
        default=os.environ.get("GITHUB_SHA"),
        help="declared source commit (not proof of binary/source matching)",
    )
    parser.add_argument(
        "--previous-report", type=Path, help="previous board JSON for progress diff"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--strict-links",
        action="store_true",
        help="fail only on missing CTest anchor names",
    )
    parser.add_argument(
        "--require-evidence",
        action="store_true",
        help="fail if any registered case lacks independent physical evidence",
    )
    args = parser.parse_args(argv)
    try:
        registry, mapping = read_json(args.registry), read_json(args.links)
        registered = read_manifest(args.ctest_manifest)
        junit = read_junits(args.junit)
        previous = (
            read_json(args.previous_report)
            if args.previous_report is not None
            else None
        )
        if previous is not None and previous.get("schema") != SCHEMA:
            raise PortfolioError("incompatible previous portfolio report")
        report, physical = make_portfolio(
            registry,
            mapping,
            registered,
            junit,
            args.evidence_dir,
            args.source_archive,
            args.source_revision,
        )
        out = args.output_dir
        md = markdown_report(report, previous)
        atomic_write(
            out / "qualification_portfolio.json",
            json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        )
        atomic_write(
            out / "physical_evidence_gate.json",
            json.dumps(physical, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        )
        atomic_write(out / "qualification_portfolio.md", md)
        write_csv(out / "qualification_portfolio.csv", report["cases"])
        summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
        if summary_path:
            with Path(summary_path).open("a", encoding="utf-8") as stream:
                stream.write(md.split("## Tests réellement fournis")[0] + "\n")
        print(
            f"CFDX_PORTFOLIO: {report['summary']['registered_cases']} registered, "
            f"{report['summary']['unregistered_backlog']} backlog, "
            f"{report['summary']['physical_evidence_gates_pass']} physical evidence passes, "
            f"{len(report['summary']['registered_missing_ctest_anchors'])} missing CTest anchors; "
            f"report={out / 'qualification_portfolio.md'}"
        )
        if args.strict_links and report["summary"]["registered_missing_ctest_anchors"]:
            return 2
        if args.require_evidence and not physical["all_evidence_gates_pass"]:
            return 2
        return 0
    except (OSError, ValueError, TypeError, KeyError, ET.ParseError) as exc:
        print(f"CFDX_PORTFOLIO: INPUT_ERROR {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
