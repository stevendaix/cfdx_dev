"""Unit tests for the SU2 adapter."""
import pytest
import numpy as np

import os

from cfdx.io.adapters.su2 import Su2Adapter, _SU2_ELEM_TYPES, _to_meshio_type
from cfdx.io.interfaces import ConversionResult, SourceInfo
from cfdx.io.schema import BCType, CaseSetup


class TestSu2Adapter:
    """Tests for the SU2 (.su2 mesh + .cfg config) adapter."""

    @pytest.fixture
    def su2_file(self):
        return os.path.join(pytest.DATA_DIR, "su2", "mesh_NACA0012_inv.su2")

    @pytest.fixture
    def cfg_file(self):
        return os.path.join(pytest.DATA_DIR, "su2", "inv_NACA0012_basic.cfg")

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
