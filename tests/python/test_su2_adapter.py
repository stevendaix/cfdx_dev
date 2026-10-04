"""Unit tests for the SU2 adapter."""
import pytest
import numpy as np
import os
from pathlib import Path

from cfdx.io.adapters.su2 import Su2Adapter, _SU2_ELEM_TYPES, _to_meshio_type
from cfdx.io.interfaces import ConversionResult, SourceInfo
from cfdx.io.schema import BCType, CaseSetup, Severity


class TestSu2Adapter:
    """Tests for the SU2 (.su2 mesh + .cfg config) adapter."""

    @pytest.fixture
    def su2_file(self):
        return os.path.join(pytest.DATA_DIR, "su2", "mesh_NACA0012_inv.su2")

    @pytest.fixture
    def cfg_file(self):
        return os.path.join(pytest.DATA_DIR, "su2", "inv_NACA0012_basic.cfg")

    @pytest.fixture
    def solution_file(self):
        return os.path.join(pytest.DATA_DIR, "su2", "mismatched_solution.csv")

    def test_adapter_metadata(self):
        adapter = Su2Adapter()
        assert adapter.solver_name == "SU2"
        assert adapter.format_name == "su2"

    def test_detect_source(self, su2_file):
        adapter = Su2Adapter()
        info = SourceInfo()
        detected = adapter.detect_source(su2_file, info)
        assert detected
        assert info.solver == "SU2"
        assert info.format == "su2"

    def test_parse_mesh(self, su2_file):
        adapter = Su2Adapter()
        ok = adapter.parse_mesh(su2_file)
        assert ok
        assert len(adapter._points) == 5233
        assert len(adapter._elements) == 10216
        assert len(adapter._boundaries) == 2
        boundary_names = [b["name"] for b in adapter._boundaries]
        assert "airfoil" in boundary_names
        assert "farfield" in boundary_names

    def test_parse_config(self, cfg_file):
        adapter = Su2Adapter()
        ok = adapter.parse_config(cfg_file)
        assert ok
        assert "MACH_NUMBER" in adapter._cfg_params
        assert float(adapter._cfg_params["MACH_NUMBER"]) == 0.8
        assert "AOA" in adapter._cfg_params
        assert adapter.setup.initial_condition.velocity > 0
        assert adapter.setup.turbulence_model == "laminar"
        assert adapter.setup.numerics.linear_solver == "fgmres"


    def test_parse_incompressible_naca0012_config(self):
        """INC_RANS must not be classified as compressible RANS."""
        cfg_file = os.path.join(pytest.DATA_DIR, "su2", "incomp_NACA0012.cfg")
        adapter = Su2Adapter()
        assert adapter.parse_config(cfg_file)

        assert adapter.setup.physics_model == "incompressible_rans"
        assert adapter.setup.energy_model == "isothermal"
        assert adapter.setup.turbulence_model == "spalart_allmaras"
        np.testing.assert_allclose(
            adapter.setup.initial_condition.velocity_vector,
            [51.36481493540834, 9.0570027322096198, 0.0],
        )
        np.testing.assert_allclose(
            adapter.setup.initial_condition.velocity,
            np.hypot(51.36481493540834, 9.0570027322096198),
        )
        assert adapter.setup.ref_density == 1.0
        assert adapter.setup.materials[0].density == 2.13163
        assert adapter.setup.materials[0].dynamic_viscosity == 1.853e-05

    def test_convert_full(self, su2_file, cfg_file):
        adapter = Su2Adapter()
        result = ConversionResult(source=SourceInfo())
        ok = adapter.convert(su2_file, result)
        assert ok
        assert not result.gap_report.has_blocking()
        assert result.mesh is not None
        assert result.mesh["points"].shape == (5233, 3)
        assert "triangle" in result.mesh["cells"]
        assert len(result.setup.boundary_conditions) == 2
        assert result.setup.mesh_info.n_vertices == 5233

    def test_convert_mesh_only(self, su2_file):
        """Test conversion when only .su2 is provided (no .cfg)."""
        adapter = Su2Adapter()
        result = ConversionResult(source=SourceInfo())
        ok = adapter.convert(su2_file, result)
        assert ok
        # Config is optional — should be non-blocking
        nonblock = result.gap_report.n_unsupported_nonblocking()
        assert nonblock >= 0

    def test_element_type_lookup(self):
        assert _SU2_ELEM_TYPES[5] == ("TRIANGLE", 3)
        assert _SU2_ELEM_TYPES[10] == ("HEXAHEDRON", 8)
        assert _SU2_ELEM_TYPES[9] == ("TETRAHEDRON", 4)

    def test_to_meshio_type(self):
        assert _to_meshio_type("TRIANGLE") == "triangle"
        assert _to_meshio_type("HEXAHEDRON") == "hexahedron"
        assert _to_meshio_type("TETRAHEDRON") == "tetra"
        assert _to_meshio_type("PRISM") == "wedge"

    def test_gap_report_has_supported_entries(self, su2_file):
        adapter = Su2Adapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(su2_file, result)
        supported = result.gap_report.n_supported()
        assert supported > 0

    def test_gap_report_documents_unsupported(self, su2_file):
        adapter = Su2Adapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(su2_file, result)
        findings = result.gap_report.findings
        categories = [f.feature for f in findings]
        assert "multiphase" in categories or "turbomachinery" in categories

    def test_parse_solution_csv(self, solution_file):
        """Test parsing SU2 solution.csv file."""
        adapter = Su2Adapter()
        scalar_fields, vec_fields = adapter.parse_solution_csv(solution_file)

        # Check scalar fields
        assert "Density" in scalar_fields
        assert "Pressure" in scalar_fields
        assert "Temperature" in scalar_fields
        assert scalar_fields["Density"].shape == (10,)
        assert scalar_fields["Pressure"].shape == (10,)
        assert scalar_fields["Temperature"].shape == (10,)
        np.testing.assert_allclose(scalar_fields["Density"], [1.225, 1.200, 1.180, 1.160, 1.140, 1.120, 1.100, 1.080, 1.060, 1.040])
        np.testing.assert_allclose(scalar_fields["Pressure"], [101325, 100000, 99000, 98000, 97000, 96000, 95000, 94000, 93000, 92000])

        # Check vector fields
        assert "Velocity" in vec_fields
        assert vec_fields["Velocity"].shape == (10, 3)
        np.testing.assert_allclose(vec_fields["Velocity"][:, 0], [50.0, 48.0, 46.0, 45.0, 44.0, 43.0, 42.0, 41.0, 40.0, 39.0])
        np.testing.assert_allclose(vec_fields["Velocity"][:, 1], [0.0, 1.0, 2.0, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5])
        np.testing.assert_allclose(vec_fields["Velocity"][:, 2], [0.0] * 10)

    def test_import_results(self, su2_file, solution_file):
        """Test importing solution.csv into ConversionResult."""
        adapter = Su2Adapter()
        result = ConversionResult(source=SourceInfo())

        # First convert mesh to populate mesh_info.n_cells
        ok = adapter.convert(su2_file, result)
        assert ok

        # Clear any fields that were imported during convert (since convert looks for solution.csv)
        result.scalar_fields.clear()
        result.vec_fields.clear()

        # Now import results explicitly
        ok = adapter.import_results(solution_file, result)
        assert not ok
        assert result.gap_report.has_blocking()
        assert result.scalar_fields == []
        assert result.vec_fields == []
        mismatch = [f for f in result.gap_report.findings if f.feature == "field_size_mismatch"]
        assert mismatch

    def test_import_results_missing_file(self, su2_file):
        """Test import_results with missing solution file."""
        adapter = Su2Adapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(su2_file, result)

        # Clear gap report to isolate the import_results call
        result.gap_report = result.gap_report.__class__()

        ok = adapter.import_results("/nonexistent/solution.csv", result)
        assert not ok
        # Should have a non-blocking gap for missing solution
        findings = result.gap_report.findings
        solution_findings = [f for f in findings if f.feature == "solution_csv"]
        assert len(solution_findings) > 0
        assert solution_findings[0].severity == Severity.UNSUPPORTED_NONBLOCK

    def test_convert_with_solution(self, su2_file, cfg_file, solution_file):
        """Test full convert including solution import."""
        # Copy solution.csv to expected location next to mesh file
        import shutil
        target = Path(su2_file).with_name(Path(su2_file).stem + "_solution.csv")
        shutil.copy(solution_file, target)

        try:
            adapter = Su2Adapter()
            result = ConversionResult(source=SourceInfo())
            ok = adapter.convert(su2_file, result)

            assert not ok
            assert result.gap_report.has_blocking()
            assert result.scalar_fields == []
            assert result.vec_fields == []
        finally:
            if target.exists():
                target.unlink()
