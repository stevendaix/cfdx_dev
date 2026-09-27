"""Unit tests for the Code_Saturne adapter."""
import pytest
import os

import numpy as np
import meshio

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
        vtk_file = tmp_path / "results.vtk"
        meshio.Mesh(
            points=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            cells=[("triangle", np.array([[0, 1, 2]]))],
            point_data={"pressure": np.array([1.0, 2.0, 3.0])},
        ).write(str(vtk_file))
        ok = adapter.import_results(str(vtk_file), result)
        assert ok
        assert not result.gap_report.has_blocking()
        assert [name for name, _ in result.scalar_fields] == ["pressure"]
        assert result.field_data_source == str(vtk_file)

    def test_import_results_unreadable_file_returns_false(self, xml_file, tmp_path):
        """A detected-but-unreadable export must NOT be reported as imported."""
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        vtk_file = tmp_path / "results.vtk"
        vtk_file.write_text("# vtk")
        ok = adapter.import_results(str(vtk_file), result)
        assert not ok
        assert not result.scalar_fields and not result.vec_fields
        assert any(
            f.feature == "neutral_vtk" and f.severity != "supported"
            for f in result.gap_report.findings
        )

    def test_import_results_detect_only_format_returns_false(self, xml_file, tmp_path):
        """MED/CGNS are detected but never parsed, so no True is reported."""
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        med_file = tmp_path / "results.med"
        med_file.write_text("dummy")
        ok = adapter.import_results(str(med_file), result)
        assert not ok
        assert any(
            f.feature == "neutral_med" and f.severity != "supported"
            for f in result.gap_report.findings
        )

    def test_import_results_missing_file_returns_false(self, xml_file, tmp_path):
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        ok = adapter.import_results(str(tmp_path / "absent.vtk"), result)
        assert not ok

    def test_python_partial_extraction_reported(self, tmp_path):
        """Non-literal case.set() calls must downgrade the status to partial."""
        py_file = tmp_path / "case.py"
        py_file.write_text(
            "from user_usr import *\n"
            'ref = "k-epsilon"\n'
            'case.set("turbulence_model", ref)\n'
            "for p in walls:\n"
            '    case.set(p, "wall")\n'
        )
        adapter = SaturneAdapter()
        adapter.parse_python_setup(str(py_file))
        assert adapter._python_extraction["status"] == "partial"
        assert adapter._python_extraction["entries"] == 0
        assert adapter._python_extraction["unresolved_case_set_calls"] == 2
        # No literal entry means no physics is invented
        assert adapter.setup.turbulence_model != "k_epsilon"

    def test_python_literal_only_extraction(self, tmp_path):
        py_file = tmp_path / "case.py"
        py_file.write_text('case.set("turbulence_model", "k-epsilon")\n')
        adapter = SaturneAdapter()
        adapter.parse_python_setup(str(py_file))
        assert adapter._python_extraction["status"] == "literal_only"
        assert adapter._python_extraction["entries"] == 1

    def test_gap_report_marks_python_as_approximated_not_supported(self, tmp_path):
        py_file = tmp_path / "case.py"
        py_file.write_text('case.set("turbulence_model", "k-epsilon")\n')
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.parse_python_setup(str(py_file))
        adapter._populate_gap_report(result.gap_report)
        py_findings = [
            f for f in result.gap_report.findings if f.feature == "saturne_py"
        ]
        assert py_findings
        assert all(f.severity != "supported" for f in py_findings)

    def test_import_results_unsupported_format(self, xml_file, tmp_path):
        adapter = SaturneAdapter()
        result = ConversionResult(source=SourceInfo())
        dat_file = tmp_path / "results.dat"
        dat_file.write_text("# some data")
        ok = adapter.import_results(str(dat_file), result)
        assert not ok
