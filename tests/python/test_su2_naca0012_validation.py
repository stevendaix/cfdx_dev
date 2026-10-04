"""Validation test for importing an official SU2 NACA0012 case.

The mesh/config pair is already present in tests/data/su2/.  The test
reconstructs a normal SU2 case directory (mesh.su2 + config.cfg) and runs
the same directory-level conversion path used by the unified importer.
It deliberately does not treat the synthetic solution.csv fixture as
validation evidence.
"""

from pathlib import Path
import shutil

from cfdx.io.adapters.su2 import Su2Adapter
from cfdx.io.interfaces import ConversionResult, SourceInfo


def test_read_su2_naca0012_validation_case(tmp_path, request):
    data_dir = Path(request.config.rootdir) / "tests" / "data" / "su2"

    mesh = data_dir / "mesh_NACA0012_inv.su2"
    cfg = data_dir / "inv_NACA0012_basic.cfg"

    assert mesh.is_file()
    assert cfg.is_file()

    case_dir = tmp_path / "naca0012"
    case_dir.mkdir()
    shutil.copy2(mesh, case_dir / "mesh.su2")
    shutil.copy2(cfg, case_dir / "config.cfg")

    adapter = Su2Adapter()
    result = ConversionResult(source=SourceInfo())

    ok = adapter.convert(case_dir, result)

    assert ok, [f"{f.feature}: {f.message}" for f in result.gap_report.findings if f.severity.value.endswith("blocking")]

    # Mesh contract from the real SU2 NACA0012 validation mesh.
    assert result.mesh is not None
    assert result.mesh["points"].shape == (5233, 3)
    assert len(result.mesh["cells"]["triangle"]) == 10216
    assert result.setup.mesh_info.dimension == 2
    assert result.setup.mesh_info.n_vertices == 5233
    assert result.setup.mesh_info.n_cells == 10216
    assert result.setup.mesh_info.n_patches == 2
    assert result.setup.mesh_info.n_boundary_faces == 250

    # Configuration contract from the SU2 validation case.
    assert result.setup.physics_model == "compressible_euler"
    assert result.setup.energy_model == "compressible"
    assert result.setup.mesh_info.n_vertices == 5233
    assert result.setup.initial_condition.pressure == 101325.0
    assert result.setup.initial_condition.temperature == 288.15
    assert result.setup.numerics.max_iterations == 1000

    assert {b.patch_name for b in result.setup.boundary_conditions} == {
        "airfoil",
        "farfield",
    }

    # No result file is supplied in this validation reader test.  The
    # adapter must report that explicitly as non-blocking rather than
    # accepting the repository's intentionally incompatible 10-row fixture.
    assert any(
        f.feature == "solution_csv" and f.severity.value.endswith("nonblocking")
        for f in result.gap_report.findings
    )
    assert not result.gap_report.has_blocking()
