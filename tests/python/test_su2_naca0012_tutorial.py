"""The public NACA reader command cannot claim a CFDX flow validation."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("naca_reader_cli", ROOT / "validation/run_su2_naca0012.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_real_fixture_reports_reader_and_production_separately(tmp_path):
    report = module.qualify_reader(tmp_path)
    assert report["reader_passed"]
    assert all(report["reader_checks"].values())
    assert not report["canonical_conversion_passed"]
    blocking = [item for item in report["canonical_conversion_findings"]["findings"]
                if item["severity"] == "unsupported_blocking"]
    assert any(item["feature"] == "canonical_selection" and "JST" in item["detail"]
               for item in blocking)
    assert report["boundary_initial_condition_equivalence_passed"] == all(report["physical_metadata_checks"].values())
    assert not report["production_supported"]
    assert report["physical_validation"] == "NOT_RUN"
    assert not list(tmp_path.glob("*.cfdx.h5"))


def test_default_command_cannot_exit_success_for_an_unrun_flow(tmp_path):
    assert module.main(["--output-dir", str(tmp_path)]) == 2
    assert module.main(["--output-dir", str(tmp_path), "--reader-only"]) == 0
