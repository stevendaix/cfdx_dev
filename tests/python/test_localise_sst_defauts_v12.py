"""Validation de la localisation SST sur vrais index de cellules du maillage."""
from __future__ import annotations

from pathlib import Path

import pytest

from scripts.generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
from scripts.localise_sst_defauts_v12 import localisation


def test_retrouve_cellule_defaut_pres_paroi(tmp_path: Path):
    mesh_path = tmp_path / "mesh.h5"
    write_case(mesh_path, make_structured_mesh("smoke"), mesh_source="fixture", max_iterations=1)
    log = tmp_path / "solver.log"
    log.write_text(
        "SST_TRANSPORT_DIAGNOSTIC iteration=100 k_local_relative=2.248040593e-01 "
        "omega_local_relative=6.073055520e-02 "
        "k_local_worst_cell=199 omega_local_worst_cell=329\n"
        "SST_TRANSPORT_DIAGNOSTIC iteration=1 k_local_relative=5 "
        "omega_local_relative=5 k_local_worst_cell=2 "
        "omega_local_worst_cell=3\n", encoding="utf-8")
    report = localisation(mesh_path, log)
    assert report["status"] == "DIAGNOSTIC_ONLY_NOT_QUALIFIED"
    assert report["iteration"] == 100
    assert report["fields"]["k"]["worst_cell"] == 199
    assert report["fields"]["k"]["defect_local_relative"] == pytest.approx(0.2248040593)
    assert report["fields"]["k"]["point_over_h"][1] > 8.9
    assert report["fields"]["omega"]["worst_cell"] == 329
    assert report["fields"]["omega"]["point_over_h"][0] == pytest.approx(1.0)
    assert report["fields"]["omega"]["defect_local_relative"] == pytest.approx(0.0607305552)


def test_rejette_labsence_de_certificat_local(tmp_path: Path):
    mesh_path = tmp_path / "mesh.h5"
    write_case(mesh_path, make_structured_mesh("smoke"), mesh_source="fixture", max_iterations=1)
    log = tmp_path / "solver.log"
    log.write_text("SST update rejected: no local residual\n")
    with pytest.raises(ValueError, match="missing SST_TRANSPORT_DIAGNOSTIC"):
        localisation(mesh_path, log)


def test_v14_egalite_iterations_choisit_la_derniere_tentative(tmp_path: Path):
    mesh_path = tmp_path / "mesh.h5"
    write_case(mesh_path, make_structured_mesh("smoke"),
               mesh_source="fixture", max_iterations=1)
    log = tmp_path / "solver.log"
    log.write_text(
        "SST_TRANSPORT_DIAGNOSTIC iteration=800 k_local_relative=2e-8 "
        "omega_local_relative=4e-9 k_local_worst_cell=210 "
        "omega_local_worst_cell=222\n"
        "SST_TRANSPORT_DIAGNOSTIC iteration=1 k_local_relative=0.4 "
        "omega_local_relative=0.2 k_local_worst_cell=3 "
        "omega_local_worst_cell=4\n"
        "SST_TRANSPORT_DIAGNOSTIC iteration=800 k_local_relative=1.5e-8 "
        "omega_local_relative=2e-9 k_local_worst_cell=211 "
        "omega_local_worst_cell=223\n", encoding="utf-8")
    sample = localisation(mesh_path, log)
    assert sample["iteration"] == 800
    assert sample["fields"]["k"]["worst_cell"] == 211
    assert sample["fields"]["omega"]["worst_cell"] == 223
    assert sample["fields"]["k"]["defect_local_relative"] == pytest.approx(1.5e-8)
