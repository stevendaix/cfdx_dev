"""Unit tests for the Fluent .cas/.dat adapter."""
import pytest
import numpy as np

import os

from cfdx.io.adapters.fluent import FluentAdapter, _parse_s_expressions, _try_pyfluent
from cfdx.io.interfaces import ConversionResult, SourceInfo
from cfdx.io.schema import BCType

PYFLUENT_AVAILABLE = _try_pyfluent()


class TestFluentAdapter:
    """Tests for the Fluent legacy ASCII .cas/.dat adapter."""

    @pytest.fixture
    def cas_file(self):
        return os.path.join(pytest.DATA_DIR, "fluent", "cavity.cas")

    def test_adapter_metadata(self):
        adapter = FluentAdapter()
        assert adapter.solver_name == "Fluent"
        assert adapter.format_name == "cas_dat"

    def test_detect_source(self, cas_file):
        adapter = FluentAdapter()
        info = SourceInfo()
        detected = adapter.detect_source(cas_file, info)
        assert detected
        assert info.solver == "Fluent"
        assert info.format == "cas_dat"

    def test_parse_cas(self, cas_file):
        adapter = FluentAdapter()
        ok = adapter.parse_cas(cas_file)
        assert ok
        assert len(adapter._points) == 30
        assert len(adapter._face_nodes) > 0
        assert len(adapter._cell_nodes) >= 1
        assert len(adapter._zones) >= 8
        assert len(adapter._materials) == 3

    def test_parse_cas_materials(self, cas_file):
        adapter = FluentAdapter()
        adapter.parse_cas(cas_file)
        mat_names = [m["name"] for m in adapter._materials]
        assert "air" in mat_names
        assert "water" in mat_names
        assert "steel" in mat_names

        air = [m for m in adapter._materials if m["name"] == "air"][0]
        assert air["density"] == 1.225

    def test_parse_cas_boundary_zones(self, cas_file):
        adapter = FluentAdapter()
        adapter.parse_cas(cas_file)
        zone_names = [z["name"] for z in adapter._zones]
        assert "wall" in zone_names
        assert "inlet" in zone_names
        assert "outlet" in zone_names
        assert "symmetry" in zone_names

    def test_convert_full(self, cas_file):
        adapter = FluentAdapter()
        result = ConversionResult(source=SourceInfo())
        ok = adapter.convert(cas_file, result)
        assert ok
        assert not result.gap_report.has_blocking()
        assert result.mesh is not None
        assert result.mesh["points"].shape == (30, 3)
        assert len(result.setup.boundary_conditions) >= 8
        assert len(result.setup.materials) == 3

    def test_boundary_conditions_mapped(self, cas_file):
        adapter = FluentAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(cas_file, result)
        bc_names = {bc.patch_name: bc.type for bc in result.setup.boundary_conditions}
        assert bc_names.get("wall") == BCType.WALL
        assert bc_names.get("inlet") == BCType.INLET
        assert bc_names.get("outlet") == BCType.PRESSURE_OUTLET

    def test_model_settings_parsed(self, cas_file):
        adapter = FluentAdapter()
        adapter.parse_cas(cas_file)
        assert len(adapter._model_settings) > 0

    def test_s_expression_parser(self):
        text = "(1 \"COMPILED\") (2 3 6\\n 0.0 0.0 0.0\\n )"
        sections = _parse_s_expressions(text)
        assert len(sections) == 2
        assert sections[0][0] == 1
        assert sections[1][0] == 2

    def test_gap_report_supported_findings(self, cas_file):
        adapter = FluentAdapter()
        result = ConversionResult(source=SourceInfo())
        adapter.convert(cas_file, result)
        assert result.gap_report.n_supported() > 0

    @pytest.mark.skipif(not PYFLUENT_AVAILABLE, reason="pyfluent not installed")
    def test_detect_cas_h5(self):
        adapter = FluentAdapter()
        info = SourceInfo()
        h5_file = os.path.join(pytest.DATA_DIR, "fluent", "h5", "mixing_elbow.cas.h5")
        detected = adapter.detect_source(h5_file, info)
        assert detected
        assert info.solver == "Fluent"
        assert info.format == "cas_h5"

    @pytest.mark.skipif(not PYFLUENT_AVAILABLE, reason="pyfluent not installed")
    def test_parse_cas_h5(self):
        adapter = FluentAdapter()
        h5_file = os.path.join(pytest.DATA_DIR, "fluent", "h5", "mixing_elbow.cas.h5")
        ok = adapter.parse_cas_h5(h5_file)
        assert ok
        assert len(adapter._points) > 0
        assert len(adapter._cell_nodes) > 0
        assert len(adapter._zones) > 1  # at least wall-inlet + wall-outlet + interior
        assert adapter.setup.source.version == "hdf5"

    @pytest.mark.skipif(not PYFLUENT_AVAILABLE, reason="pyfluent not installed")
    def test_convert_cas_h5_full(self):
        adapter = FluentAdapter()
        h5_file = os.path.join(pytest.DATA_DIR, "fluent", "h5", "mixing_elbow.cas.h5")
        result = ConversionResult(source=SourceInfo())
        ok = adapter.convert(h5_file, result)
        assert ok
        assert not result.gap_report.has_blocking()
        assert result.mesh is not None
        assert result.mesh["points"].shape[1] == 3  # 3D mesh
        assert result.setup.mesh_info.n_vertices > 0
        # Check materials are extracted
        assert len(result.setup.materials) >= 1
        # Check turbulence model is extracted
        assert result.setup.turbulence_model  # non-empty
        # Check reference values
        assert result.setup.ref_density == 1.225
        assert result.setup.ref_velocity == 1.0
        # Check boundary conditions are created from zones
        assert len(result.setup.boundary_conditions) >= 1
        # Check energy model
        assert result.setup.energy_model == "isothermal"

    @pytest.mark.skipif(not PYFLUENT_AVAILABLE, reason="pyfluent not installed")
    def test_import_dat_h5(self):
        adapter = FluentAdapter()
        h5_file = os.path.join(pytest.DATA_DIR, "fluent", "h5", "mixing_elbow.cas.h5")
        dat_file = os.path.join(pytest.DATA_DIR, "fluent", "h5", "mixing_elbow.dat.h5")
        if not os.path.exists(dat_file):
            pytest.skip("mixing_elbow.dat.h5 not available")
        adapter.parse_cas_h5(h5_file)
        result = ConversionResult(source=SourceInfo())
        ok = adapter.import_results(dat_file, result)
        assert ok
        assert "scalar" in result.available_fields
        assert "vector" in result.available_fields
        assert len(result.available_fields["scalar"]) >= 1
