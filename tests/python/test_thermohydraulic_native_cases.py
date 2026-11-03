"""Native CFDX TH campaign and wall diagnostic contracts (no fake physical PASS)."""

from __future__ import annotations

import copy
import csv
import json
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts import run_thermohydraulic_cases as campaign  # noqa: E402
from scripts import extract_th03_wall_diagnostics as wall  # noqa: E402


def test_example_campaign_reports_all_unavailable_critical_physics() -> None:
    path=ROOT/"benchmarks/thermohydraulics/campaign.example.json"
    data=json.loads(path.read_text(encoding="utf-8"))
    report=campaign.preflight(data,path.parent,None)
    assert report["status"] == "BLOCKED"
    assert [r["id"] for r in report["cases"]] == list(campaign.SCENARIOS)
    by_id={c["id"]:c for c in report["cases"]}
    assert any("periodic" in x for x in by_id["TH-01"]["blockers"])
    assert any("SST+energy+buoyancy" in x for x in by_id["TH-07"]["blockers"])
    assert any("transient" in x for x in by_id["TH-08"]["blockers"])
    assert all(r["physical_qualification"] == "NOT_QUALIFIED" for r in report["cases"])


def test_three_mesh_scenario_is_mandatory(tmp_path: Path) -> None:
    data={"schema_version":1,"cases":[{"id":"TH-03","model":"RANS_SST",
          "scenarios":[{"id":"Re_h_37500","runs":[{"mesh_id":"fine"}]}]}]}
    with pytest.raises(campaign.CampaignError,match="three mesh"):
        campaign.preflight(data,tmp_path,None)


def test_reject_path_escape_and_wrong_mesh_format(tmp_path: Path) -> None:
    path=ROOT/"benchmarks/thermohydraulics/campaign.example.json"
    data=json.loads(path.read_text(encoding="utf-8"))
    data["cases"]=[copy.deepcopy(data["cases"][1])]
    first=data["cases"][0]["scenarios"][0]["runs"][0]
    first["case_file"]="../../../etc/passwd.cfdx.h5"
    with pytest.raises(campaign.CampaignError,match="escapes"):
        campaign.preflight(data,tmp_path,None)
    first["case_file"]="mesh.msh"
    with pytest.raises(campaign.CampaignError,match=".cfdx.h5"):
        campaign.preflight(data,tmp_path,None)


def test_wall_reattachment_signed_zero_crossing(tmp_path: Path) -> None:
    path=tmp_path/"face_fluxes.csv"
    with path.open("w",newline="",encoding="utf-8") as target:
        writer=csv.DictWriter(target,fieldnames=["patch","x_m","area_m2",
            "wall_shear_x_Pa","wall_shear_method"])
        writer.writeheader()
        # negative at x/h=2, positive at x/h=4 -> exactly 3
        for x,tau in [(1,-3),(2,-2),(4,2),(6,4)]:
            writer.writerow({"patch":"lower_wall","x_m":x*0.01,
                "area_m2":0.2,"wall_shear_x_Pa":tau,
                "wall_shear_method":"first_order_one_sided_diagnostic"})
    profiles,root=wall.calculate_wall_diagnostics(path,patch="lower_wall",
        density=1,bulk_speed=4,step_height=0.01)
    assert len(profiles)==4
    assert root == pytest.approx(3)
    assert profiles[0]["value"] < 0 < profiles[-1]["value"]


def test_wall_diagnostic_fails_when_no_reattachment_or_bad_method(tmp_path: Path) -> None:
    path=tmp_path/"face_fluxes.csv"
    path.write_text("patch,x_m,area_m2,wall_shear_x_Pa,wall_shear_method\n"
        "lower_wall,0.01,1,-3,first_order_one_sided_diagnostic\n"
        "lower_wall,0.02,1,-1,first_order_one_sided_diagnostic\n"
        "lower_wall,0.03,1,-0.2,first_order_one_sided_diagnostic\n",encoding="utf-8")
    with pytest.raises(wall.WallDiagnosticError,match="crossing"):
        wall.calculate_wall_diagnostics(path,patch="lower_wall",
            density=1,bulk_speed=1,step_height=0.01)
    content=path.read_text(encoding="utf-8").replace(
        "first_order_one_sided_diagnostic","unqualified_mystery_method")
    path.write_text(content,encoding="utf-8")
    with pytest.raises(wall.WallDiagnosticError,match="reconstruction missing"):
        wall.calculate_wall_diagnostics(path,patch="lower_wall",
            density=1,bulk_speed=1,step_height=0.01)
