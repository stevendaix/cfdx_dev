"""The A/B diagnostic must never be mistaken for a qualification."""
import csv
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/naca_initialization_ab.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("naca_initialization_ab", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


def create_history(path: Path, momentum: float, scheme: str = "TVD_VANALBADA"):
    from naca_convergence_diagnosis import REQUIRED
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=REQUIRED)
        writer.writeheader()
        for i in range(1, 3):
            row = dict.fromkeys(REQUIRED, "0")
            row.update(
                iteration=str(i),
                algorithm="SIMPLE",
                convection_scheme=scheme,
                high_order_fraction="1",
                momentum_rhs_relative=str(momentum / i),
                momentum_x_rhs_relative=str(momentum / i),
                momentum_y_rhs_relative="0.01",
                momentum_component_max="0.01",
            )
            writer.writerow(row)


def test_valid_matching_mesh_history_is_diagnostic_only(tmp_path):
    cold, warm = tmp_path / "cold.csv", tmp_path / "warm.csv"
    create_history(cold, 0.01)
    create_history(warm, 0.003)
    mesh = tmp_path / "mesh.h5"
    mesh.write_bytes(b"mesh-content")
    report = module.compare(cold, warm, mesh, mesh)
    assert report["common_iterations"] == 2
    assert report["seed_independence_proven"] is False
    assert (
        report["metrics"]["momentum_rhs_relative"]["cold_over_warm"]
        == pytest.approx(10 / 3)
    )


def test_reject_different_mesh_bytes(tmp_path):
    cold, warm = tmp_path / "cold.csv", tmp_path / "warm.csv"
    create_history(cold, 0.01)
    create_history(warm, 0.003)
    mesh_c, mesh_w = tmp_path / "mesh-c", tmp_path / "mesh-w"
    mesh_c.write_bytes(b"a")
    mesh_w.write_bytes(b"b")
    with pytest.raises(ValueError, match="identical mesh bytes"):
        module.compare(cold, warm, mesh_c, mesh_w)


def test_reject_different_schemes(tmp_path):
    cold, warm = tmp_path / "cold.csv", tmp_path / "warm.csv"
    create_history(cold, 0.01)
    create_history(warm, 0.003, "UPWIND")
    mesh = tmp_path / "mesh"
    mesh.write_bytes(b"a")
    with pytest.raises(ValueError, match="different convection_scheme"):
        module.compare(cold, warm, mesh, mesh)


def test_reject_negative_residual(tmp_path):
    cold, warm = tmp_path / "cold.csv", tmp_path / "warm.csv"
    create_history(cold, -0.01)
    create_history(warm, 0.003)
    mesh = tmp_path / "mesh"
    mesh.write_bytes(b"a")
    with pytest.raises(ValueError, match="invalid nonnegative metric"):
        module.compare(cold, warm, mesh, mesh)