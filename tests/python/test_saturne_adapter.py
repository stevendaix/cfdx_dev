"""Unit tests for the Code_Saturne adapter."""
import pytest
import os

import os

from cfdx.io.adapters.saturne import SaturneAdapter
from cfdx.io.interfaces import ConversionResult, SourceInfo
from cfdx.io.schema import BCType


class TestSaturneAdapter:
    """Tests for the Code_Saturne (XML + Python) adapter."""

    @pytest.fixture
    def xml_file(self):
        return os.path.join(pytest.DATA_DIR, "saturne", "test_case.xml")

    @pytest.fixture
    def py_file(self):
        return os.path.join(pytest.DATA_DIR, "saturne", "test_case.py")

    def test_adapter_metadata(self):
        adapter = SaturneAdapter()
        assert adapter.solver_name == "Code_Saturne"
        assert adapter.format_name == "saturne_xml_py"

    def test_detect_source_xml(self, xml_file):
        adapter = SaturneAdapter()
        info = SourceInfo()
        detected = adapter.detect_source(xml_file, info)
        assert detected
        assert info.solver == "Code_Saturne"

    def test_parse_xml_setup(self, xml_file):
        adapter = SaturneAdapter()
        ok = adapter.parse_xml_setup(xml_file)
        assert ok
        assert len(adapter._entries) > 0
        assert any(e["key"] == "physics_model" for e in adapter._entries)

    def test_parse_xml_boundary_conditions(self, xml_file):
        adapter = SaturneAdapter()
        adapter.parse_xml_setup(xml_file)
        bc_names = [bc.patch_name for bc in adapter.setup.boundary_conditions]
        assert "inlet" in bc_names
        assert "outlet" in bc_names
        assert "wall" in bc_names
        assert "symmetry" in bc_names

    def test_parse_xml_materials(self, xml_file):
        adapter = SaturneAdapter()
        adapter.parse_xml_setup(xml_file)
        assert len(adapter.setup.materials) >= 1
        assert adapter.setup.materials[0].name == "air"

    def test_parse_xml_physics(self, xml_file):
        adapter = SaturneAdapter()
        adapter.parse_xml_setup(xml_file)
        assert adapter.setup.turbulence_model == "k_epsilon"
        assert adapter.setup.physics_model == "compressible_navier_stokes"

    def test_parse_python_setup(self, py_file):
        adapter = SaturneAdapter()
        ok = adapter.parse_python_setup(py_file)
        assert ok
        assert len(adapter._entries) > 0
        py_entries = [e for e in adapter._entries if e["section"] == "python"]
        assert len(py_entries) > 0

    def test_convert_blocks_without_neutral_mesh(self, xml_file):
        """Conversion must be blocked when no neutral mesh export is present."""
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        ok = adapter.convert(xml_file, result)
        assert not ok
        assert result.gap_report.has_blocking()

    def test_convert_with_synthetic_py(self, py_file):
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.parse_python_setup(py_file)
        assert adapter.setup.turbulence_model == "k_epsilon"

    def test_import_results_neutral_formats(self, xml_file, tmp_path):
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        # Create a fake .vtk file
        vtk_file = tmp_path / "results.vtk"
        vtk_file.write_text("# vtk")
        ok = adapter.import_results(str(vtk_file), result)
        assert ok
        assert not result.gap_report.has_blocking()

    def test_import_results_unsupported_format(self, xml_file, tmp_path):
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        dat_file = tmp_path / "results.dat"
        dat_file.write_text("# some data")
        ok = adapter.import_results(str(dat_file), result)
        assert not ok
