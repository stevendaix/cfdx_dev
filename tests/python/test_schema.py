"""Unit tests for the CFDX I/O schema, interfaces, and GapAnalysis."""
import pytest
import numpy as np

from cfdx.io.schema import (
    CaseSetup,
    BoundarySpec,
    MaterialSpec,
    SourceInfo,
    MeshMetadata,
    BCType,
    BCValueType,
    CFDX_SCHEMA_VERSION,
    Severity,
)
from cfdx.io.interfaces import (
    SolverAdapter,
    ConversionResult,
)
from cfdx.io.gap_analysis import GapAnalysis, GapAnalysisReport, Finding
from cfdx.io.mapping_rules import MappingRules


class TestSchema:
    def test_schema_version(self):
        assert CFDX_SCHEMA_VERSION == 1

    def test_case_setup_defaults(self):
        cs = CaseSetup()
        assert cs.schema_version == CFDX_SCHEMA_VERSION
        assert cs.mesh_info.dimension == 3
        assert cs.physics_model == "incompressible_laminar"
        assert cs.turbulence_model == "laminar"
        assert cs.energy_model == "isothermal"
        assert cs.transient is False

    def test_bc_types(self):
        assert BCType.WALL == "wall"
        assert BCType.INLET == "inlet"
        assert BCType.OUTLET == "outlet"

    def test_bc_value_types(self):
        assert BCValueType.FIXED == "fixed"
        assert BCValueType.WALL_NO_SLIP == "wall_no_slip"

    def test_find_boundary(self):
        bc1 = BoundarySpec(patch_name="inlet", type=BCType.INLET)
        bc2 = BoundarySpec(patch_name="wall", type=BCType.WALL)
        cs = CaseSetup(boundary_conditions=[bc1, bc2])
        assert cs.find_boundary("inlet") is bc1
        assert cs.find_boundary("wall") is bc2
        assert cs.find_boundary("nonexistent") is None

    def test_find_material(self):
        mat = MaterialSpec(name="air", density=1.225)
        cs = CaseSetup(materials=[mat])
        assert cs.find_material("air") is mat
        assert cs.find_material("water") is None


class TestMappingRules:
    def test_load_default(self):
        rules = MappingRules.load_default()
        assert rules.n_bc_mappings() > 10
        assert rules.n_turbulence_mappings() > 5
        assert rules.n_scheme_mappings() > 5

    def test_map_boundary_type(self):
        rules = MappingRules.load_default()
        assert rules.map_boundary_type("velocity-inlet") == BCType.INLET
        assert rules.map_boundary_type("wall") == BCType.WALL
        assert rules.map_boundary_type("symmetry") == BCType.SYMMETRY
        assert rules.map_boundary_type("unknown_type") == BCType.UNKNOWN

    def test_map_bc_value_type(self):
        rules = MappingRules.load_default()
        assert rules.map_bc_value_type("velocity-inlet") == BCValueType.FIXED
        assert rules.map_bc_value_type("wall") == BCValueType.WALL_NO_SLIP

    def test_map_turbulence(self):
        rules = MappingRules.load_default()
        assert rules.map_turbulence("NONE") == "laminar"
        assert rules.map_turbulence("SA") == "spalart_allmaras"
        assert rules.map_turbulence("KEPSILON") == "k_epsilon"

    def test_map_scheme(self):
        rules = MappingRules.load_default()
        assert rules.map_scheme("MUSCL") == "MUSCL"
        assert rules.map_scheme("GREEN_GAUSS_CELL_VOLUME") == "green_gauss_cell"

    def test_custom_mapping(self):
        rules = MappingRules()
        rules.add_bc_mapping("my_bc", BCType.INLET)
        assert rules.map_boundary_type("my_bc") == BCType.INLET


class TestGapAnalysis:
    def test_finding_creation(self):
        f = Finding(
            severity=Severity.SUPPORTED,
            category="mesh",
            feature="points",
            detail="test",
        )
        assert f.category == "mesh"
        assert f.feature == "points"

    def test_gap_analysis_no_findings(self):
        gap = GapAnalysis()
        assert len(gap.findings) == 0
        assert not gap.has_blocking()

    def test_add_supported(self):
        gap = GapAnalysis()
        gap.supported("mesh", "points", "Parsed 100 points")
        assert gap.n_supported() == 1
        assert not gap.has_blocking()

    def test_add_blocking(self):
        gap = GapAnalysis()
        gap.unsupported_blocking("mesh", "topology", "Not supported", "Use export")
        assert gap.n_unsupported_blocking() == 1
        assert gap.has_blocking()

    def test_add_nonblocking(self):
        gap = GapAnalysis()
        gap.unsupported_nonblocking("mesh", "feature", "Not supported", "Workaround")
        assert gap.n_unsupported_nonblocking() == 1
        assert not gap.has_blocking()

    def test_add_approximated(self):
        gap = GapAnalysis()
        gap.approximated("physics", "model", "Approximated", "Verify manually")
        assert gap.n_approximated() == 1
        assert not gap.has_blocking()

    def test_add_unavailable(self):
        gap = GapAnalysis()
        gap.unavailable("config", "file", "File not found")
        assert gap.n_unavailable() == 1
        assert not gap.has_blocking()

    def test_json_report(self):
        gap = GapAnalysis()
        gap.supported("mesh", "points", "Parsed")
        gap.unsupported_nonblocking("features", "x", "Not supported", "Fix it")
        report = GapAnalysisReport(gap)
        json_str = report.to_json()
        assert '"has_blocking": false' in json_str
        assert '"supported": 1' in json_str

    def test_markdown_report(self):
        gap = GapAnalysis()
        gap.supported("mesh", "points", "Parsed")
        gap.unsupported_blocking("mesh", "topology", "Blocked", "Use export")
        gap.unsupported_nonblocking("features", "x", "Not supported", "Fix")
        report = GapAnalysisReport(gap)
        md = report.to_markdown()
        assert "# Gap Analysis Report" in md
        assert "## Summary" in md
        assert "Blocking" in md


class TestConversionResult:
    def test_defaults(self):
        info = SourceInfo()
        result = ConversionResult(source=info)
        assert result.mesh is None
        assert result.setup is None
        assert len(result.scalar_fields) == 0
        assert len(result.vec_fields) == 0
        assert not result.gap_report.has_blocking()
        assert result.success
        assert not result.has_mesh
