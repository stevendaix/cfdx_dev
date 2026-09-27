"""Unit tests for the STAR-CCM+ .sim adapter."""
import pytest

import os

from cfdx.io.adapters.starccm import StarCCMAdapter, _extract_ascii_strings
from cfdx.io.interfaces import ConversionResult, SourceInfo


class TestStarCCMAdapter:
    """Tests for the STAR-CCM+ .sim adapter.

    Per issue #425: .sim binary mesh topology is NOT reconstructed.
    Tests verify that metadata extraction works and blocking gaps
    are properly documented.
    """

    @pytest.fixture
    def sim_file(self):
        return os.path.join(pytest.DATA_DIR, "starccm", "test_case.sim")

    def test_adapter_metadata(self):
        adapter = StarCCMAdapter()
        assert adapter.solver_name == "STAR-CCM+"
        assert adapter.format_name == "sim"

    def test_detect_source(self, sim_file):
        adapter = StarCCMAdapter()
        info = SourceInfo()
        detected = adapter.detect_source(sim_file, info)
        assert detected
        assert info.solver == "STAR-CCM+"
        assert info.format == "sim"

    def test_convert_blocks_on_no_binary_topology(self, sim_file):
        """Conversion must be blocked because .sim topology is not reconstructed."""
        adapter = StarCCMAdapter()
        result = ConversionResult(source=SourceInfo())
        ok = adapter.convert(sim_file, result)
        assert not ok  # blocked
        assert result.gap_report.has_blocking()

    def test_boundary_zones_detected(self, sim_file):
        adapter = StarCCMAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(sim_file, result)
        zone_names = [z["name"] for z in adapter._zones]
        assert "wall" in zone_names
        assert "inlet" in zone_names
        assert "outlet" in zone_names

    def test_materials_detected(self, sim_file):
        adapter = StarCCMAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(sim_file, result)
        mat_names = [m["name"] for m in adapter._materials]
        assert "air" in mat_names

    def test_physics_detected(self, sim_file):
        adapter = StarCCMAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(sim_file, result)
        assert "k-epsilon" in adapter._physics

    def test_import_results_unsupported(self, sim_file):
        adapter = StarCCMAdapter()
        result = ConversionResult(source=SourceInfo())
        ok = adapter.import_results(sim_file, result)
        assert not ok
        assert result.gap_report.has_blocking()

    def test_extract_ascii_strings(self):
        data = b"hello\x00world\x01foo bar baz\x00"
        strings = _extract_ascii_strings(data, min_len=3)
        assert "hello" in strings
        assert "world" in strings
        assert "foo bar baz" in strings

    def test_gap_report_documents_limitations(self, sim_file):
        adapter = StarCCMAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(sim_file, result)
        findings = result.gap_report.findings
        features = [f.feature for f in findings]
        assert "sim_topology" in features
        assert "field_data" in features
