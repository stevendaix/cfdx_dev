import csv
import json
import pathlib
import sys

import pytest

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / "scripts"))
from compare_naca_multisecant import EvidenceError, compare


def create(tmp_path, method="OFF", count=5):
    (tmp_path / "app").write_bytes(b"test binary")
    (tmp_path / "mesh").write_bytes(b"test mesh")
    fields = [
        "iteration",
        "continuity_normalized",
        "momentum_rhs_relative",
        "mass_normalized_imbalance",
        "mass_local_linf",
        "face_momentum_x_relative",
        "face_momentum_y_relative",
        "momentum_x_rhs_relative",
        "momentum_y_rhs_relative",
        "momentum_component_max",
        "velocity_change_inf",
        "pressure_change_inf",
        "pressure_residual",
        "alpha_u",
        "alpha_p",
        "high_order_fraction",
        "multisecant_relative_step",
        "multisecant_estimated_objective",
        "nonfinite_mass_faces",
        "nonfinite_velocities",
        "multisecant_proposed",
        "multisecant_accepted",
        "convection_scheme",
        "algorithm",
        "multisecant_method",
    ]
    for name, actual_method in [("base", "OFF"), ("candidate", method)]:
        with (tmp_path / f"{name}.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for i in range(1, count + 1):
                row = {key: "0.00001" for key in fields}
                row.update(
                    iteration=i,
                    alpha_u="0.6",
                    alpha_p="0.25",
                    high_order_fraction="1",
                    convection_scheme="TVD",
                    algorithm="SIMPLE",
                    multisecant_method=actual_method,
                    nonfinite_mass_faces="0",
                    nonfinite_velocities="0",
                    multisecant_proposed="0",
                    multisecant_accepted="0",
                )
                writer.writerow(row)
    manifest = {
        "binary": "app",
        "mesh": "mesh",
        "alpha_u": 0.6,
        "alpha_p": 0.25,
        "iterations": count,
        "tolerance": 1e-8,
        "runs": [
            {"name": "base", "history": "base.csv", "method": "OFF"},
            {"name": "candidate", "history": "candidate.csv", "method": method},
        ],
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    return path


def edit_row(tmp_path, key, value):
    path = tmp_path / "candidate.csv"
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
        fields = list(rows[0])
    rows[2][key] = value
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_accepts_only_as_evidence(tmp_path):
    result = compare(create(tmp_path, "ANDERSON"))
    assert result["runs"][1]["complete"]
    assert result["runs"][1]["qualification"] == "NOT_CERTIFIED"
    assert result["runs"][1]["requires_official_certificate"]


@pytest.mark.parametrize(
    "field,corruption",
    [
        ("convection_scheme", "UPWIND"),
        ("high_order_fraction", "0.5"),
        ("momentum_x_rhs_relative", "NaN"),
        ("iteration", "5"),
        ("alpha_p", "0.20"),
        ("multisecant_method", "NGMRES"),
        ("nonfinite_mass_faces", "1"),
        ("multisecant_accepted", "2"),
    ],
)
def test_bad_history_fails_closed(tmp_path, field, corruption):
    manifest = create(tmp_path, "ANDERSON")
    edit_row(tmp_path, field, corruption)
    with pytest.raises(EvidenceError):
        compare(manifest)


def test_incomplete_without_explicit_interruption_fails(tmp_path):
    manifest = create(tmp_path, "NGMRES", 3)
    data = json.loads(manifest.read_text())
    data["iterations"] = 4
    manifest.write_text(json.dumps(data))
    with pytest.raises(EvidenceError):
        compare(manifest)
    data["runs"][0]["interrupted"] = True
    data["runs"][1]["interrupted"] = True
    manifest.write_text(json.dumps(data))
    assert not compare(manifest)["runs"][1]["complete"]


def test_checksum_fails_closed(tmp_path):
    manifest = create(tmp_path, "ANDERSON")
    data = json.loads(manifest.read_text())
    data["expected_binary_sha256"] = "0" * 64
    manifest.write_text(json.dumps(data))
    with pytest.raises(EvidenceError):
        compare(manifest)


def test_relaxation_of_phys_gate_refused(tmp_path):
    manifest = create(tmp_path)
    data = json.loads(manifest.read_text())
    data["tolerance"] = 1e-5
    manifest.write_text(json.dumps(data))
    with pytest.raises(EvidenceError):
        compare(manifest)
