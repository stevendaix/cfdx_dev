"""Unit tests for the unified CFDX conversion API (cfdx.io.converter)."""

import json
import os
import shutil

import numpy as np
import pytest

from cfdx.io.converter import (
    ConversionResult,
    SolverConverter,
    available_converters,
    convert,
    detect_solver,
    get_converter,
    register_converter,
)
from cfdx.io.gap_analysis import GapAnalysis
from cfdx.io.schema import CaseSetup, Severity, SourceInfo

DATA_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "data")
)


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------


def _foam_header(obj, cls):
    return (
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
        f"    class {cls};\n    object {obj};\n}}\n"
    )


@pytest.fixture
def openfoam_case(tmp_path):
    """Minimal but valid OpenFOAM case directory (one hex cell)."""
    case = tmp_path / "foam_case"
    poly = case / "constant" / "polyMesh"
    poly.mkdir(parents=True)

    points = [
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    ]
    faces = [
        (0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
        (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
    ]
    (poly / "points").write_text(
        _foam_header("points", "vectorField")
        + f"{len(points)}\n(\n"
        + "".join(f"({' '.join(map(str, p))})\n" for p in points)
        + ")\n"
    )
    (poly / "faces").write_text(
        _foam_header("faces", "faceList")
        + f"{len(faces)}\n(\n"
        + "".join(f"{len(f)}({' '.join(map(str, f))})\n" for f in faces)
        + ")\n"
    )
    (poly / "owner").write_text(
        _foam_header("owner", "labelList") + "6\n(\n0\n1\n0\n0\n0\n0\n)\n"
    )
    (poly / "neighbour").write_text(
        _foam_header("neighbour", "labelList") + "4\n(\n1\n2\n3\n4\n)\n"
    )
    (poly / "boundary").write_text(
        _foam_header("boundary", "polyBoundaryMesh")
        + "2\n(\n"
        "    inlet\n    {\n        type            patch;\n"
        "        nFaces          1;\n        startFace       0;\n    }\n"
        "    outlet\n    {\n        type            wall;\n"
        "        nFaces          1;\n        startFace       1;\n    }\n"
        ")\n"
    )
    return case


@pytest.fixture
def su2_case(tmp_path):
    """SU2 mesh + config copied into a temp dir with matching stems.

    No solution file is staged next to the mesh: the repository fixture holds
    10 rows against 10216 cells, which is a deliberate blocking cardinality
    mismatch. Result-import handling is covered by tests/python/test_su2_adapter.py.
    """
    src = os.path.join(DATA_DIR, "su2")
    case = tmp_path / "case.su2"
    shutil.copy(os.path.join(src, "mesh_NACA0012_inv.su2"), case)
    shutil.copy(os.path.join(src, "inv_NACA0012_basic.cfg"), tmp_path / "case.cfg")
    return case


# --------------------------------------------------------------------------
# Registry / get_converter
# --------------------------------------------------------------------------


class TestRegistry:
    def test_get_converter_fluent(self):
        converter = get_converter("fluent")
        assert isinstance(converter, SolverConverter)
        assert converter.solver_name == "Fluent"
        assert ".cas" in converter.extensions
        assert ".cas.h5" in converter.extensions

    @pytest.mark.parametrize(
        "name,solver",
        [
            ("fluent", "Fluent"),
            ("openfoam", "OpenFOAM"),
            ("su2", "SU2"),
            ("starccm", "STAR-CCM+"),
            ("saturne", "Code_Saturne"),
            ("cfdx", "CFDX"),
        ],
    )
    def test_get_converter_all(self, name, solver):
        assert get_converter(name).solver_name == solver

    def test_get_converter_is_case_insensitive(self):
        assert get_converter("Fluent").solver_name == "Fluent"
        assert get_converter("OpenFOAM").solver_name == "OpenFOAM"
        assert get_converter("Code_Saturne").solver_name == "Code_Saturne"

    def test_get_converter_returns_fresh_instance(self):
        assert get_converter("su2") is not get_converter("su2")

    def test_unknown_solver_raises(self):
        with pytest.raises(ValueError):
            get_converter("not-a-solver")

    def test_available_converters_lists_all(self):
        names = available_converters()
        for expected in ("fluent", "openfoam", "su2", "starccm", "saturne"):
            assert expected in names

    def test_register_converter(self, tmp_path):
        class DummyConverter(SolverConverter):
            solver_name = "Dummy"
            format_name = "dummy"

            def _read(self, path):
                return ConversionResult(source=SourceInfo(solver="Dummy"))

        register_converter("dummy", DummyConverter)
        try:
            assert isinstance(get_converter("dummy"), DummyConverter)
        finally:
            from cfdx.io import converter as converter_module

            converter_module._CONVERTER_REGISTRY.pop("dummy", None)


# --------------------------------------------------------------------------
# Format detection
# --------------------------------------------------------------------------


class TestDetection:
    @pytest.mark.parametrize(
        "filename,expected",
        [
            ("case.cas", "fluent"),
            ("case.dat", "fluent"),
            ("case.cas.h5", "fluent"),
            ("case.dat.h5", "fluent"),
            ("case.su2", "su2"),
            ("case.cfg", "su2"),
            ("case.sim", "starccm"),
            ("case.xml", "saturne"),
            ("case.py", "saturne"),
            ("mesh.med", "saturne"),
            ("mesh.cgns", "saturne"),
            ("case.cfdx.h5", "cfdx"),
        ],
    )
    def test_extension_detection(self, filename, expected):
        assert detect_solver(filename) == expected

    def test_directory_with_polymesh_is_openfoam(self, openfoam_case):
        assert detect_solver(openfoam_case) == "openfoam"

    def test_directory_with_su2_mesh(self, tmp_path):
        (tmp_path / "mesh.su2").write_text("NDIME= 3\n")
        assert detect_solver(tmp_path) == "su2"

    def test_plain_directory_defaults_to_openfoam(self, tmp_path):
        assert detect_solver(tmp_path) == "openfoam"

    def test_unknown_extension_raises_on_convert(self, tmp_path):
        bogus = tmp_path / "mystery.zzz"
        bogus.write_text("nothing useful")
        with pytest.raises(ValueError):
            convert(bogus)

    def test_explicit_solver_overrides_detection(self, tmp_path, su2_case):
        result = convert(su2_case, solver="su2")
        assert result.source.solver == "SU2"


# --------------------------------------------------------------------------
# ConversionResult contract
# --------------------------------------------------------------------------


class TestConversionResult:
    def test_default_result_has_all_attributes(self):
        result = ConversionResult()
        assert result.source is not None
        assert result.case is None
        assert result.mesh is None
        assert result.fields == []
        assert result.boundaries == []
        assert result.materials == []
        assert result.solver_settings == {}
        assert isinstance(result.gap_analysis, GapAnalysis)
        assert result.output_path is None

    def test_success_reflects_blocking_findings(self):
        result = ConversionResult()
        assert result.success is True
        result.gap_analysis.unsupported_blocking("mesh", "topology", "no topology")
        assert result.success is False

    def test_nonblocking_findings_keep_success(self):
        result = ConversionResult()
        result.gap_analysis.unsupported_nonblocking("bc", "roughness", "dropped")
        assert result.success is True

    def test_helpers(self):
        result = ConversionResult()
        result.fields.append(("p", np.array([1.0, 2.0])))
        assert result.field_names == ["p"]
        assert np.array_equal(result.find_field("p"), [1.0, 2.0])
        assert result.find_field("missing") is None
        assert result.has_mesh is False

    def test_summary_is_json_serialisable(self):
        result = ConversionResult(source=SourceInfo(solver="Fluent"))
        json.dumps(result.summary())


# --------------------------------------------------------------------------
# convert() end to end
# --------------------------------------------------------------------------


class TestConvert:
    def test_convert_su2_autodetect(self, su2_case):
        result = convert(su2_case)
        assert result.source.solver == "SU2"
        assert result.success is True
        assert result.has_mesh is True
        assert isinstance(result.case, CaseSetup)
        assert len(result.mesh["points"]) > 0
        assert len(result.fields) > 0
        assert "airfoil" in result.boundary_names
        assert result.solver_settings["turbulence_model"]

    def test_convert_openfoam_autodetect(self, openfoam_case):
        result = convert(openfoam_case)
        assert result.source.solver == "OpenFOAM"
        assert result.success is True
        assert result.has_mesh is True
        assert set(result.boundary_names) == {"inlet", "outlet"}

    def test_convert_fluent_autodetect(self):
        result = convert(os.path.join(DATA_DIR, "fluent", "cavity.cas"))
        assert result.source.solver == "Fluent"
        assert result.has_mesh is True
        assert "air" in result.material_names
        assert "wall" in result.boundary_names

    def test_starccm_records_blocking_gap(self):
        result = convert(os.path.join(DATA_DIR, "starccm", "test_case.sim"))
        assert result.source.solver == "STAR-CCM+"
        assert result.success is False
        assert result.gap_analysis.n_unsupported_blocking() > 0

    def test_get_converter_convert_signature(self, su2_case):
        converter = get_converter("su2")
        result = converter.convert(source=su2_case)
        assert result.has_mesh is True
        assert result.success is True

    def test_missing_source_path_is_blocking(self, tmp_path):
        result = get_converter("su2").convert(tmp_path / "nope.su2")
        assert result.success is False
        assert result.gap_analysis.n_unsupported_blocking() == 1
        assert result.gap_analysis.findings[0].category == "source"

    def test_adapter_exception_is_recorded_not_raised(self, tmp_path):
        class ExplodingAdapter:
            solver_name = "Exploding"
            format_name = "boom"

            def detect_source(self, path, info):
                return False

            def convert(self, path, result):
                raise RuntimeError("boom")

        from cfdx.io.converter import AdapterConverter

        class ExplodingConverter(AdapterConverter):
            solver_name = "Exploding"
            format_name = "boom"

            def _make_adapter(self):
                return ExplodingAdapter()

        source = tmp_path / "boom.case"
        source.write_text("x")
        result = ExplodingConverter().convert(source)

        assert result.success is False
        blocking = [
            f
            for f in result.gap_analysis.findings
            if f.severity == Severity.UNSUPPORTED_BLOCK
        ]
        assert len(blocking) == 1
        assert "RuntimeError: boom" in blocking[0].detail

    def test_directory_rejected_by_file_only_converter(self, tmp_path):
        result = get_converter("starccm").read(tmp_path)
        assert result.success is False
        assert any(
            f.category == "source" for f in result.gap_analysis.findings
        )


# --------------------------------------------------------------------------
# Output writing / round trip
# --------------------------------------------------------------------------


class TestOutput:
    def test_writes_hdf5_and_reports(self, su2_case, tmp_path):
        out = tmp_path / "case.cfdx.h5"
        result = convert(su2_case, output=out)

        assert out.exists()
        assert result.output_path == str(out)
        assert (tmp_path / "case.cfdx_gap_analysis.md").exists()
        assert (tmp_path / "case.cfdx_gap_analysis.json").exists()

        report = json.loads((tmp_path / "case.cfdx_gap_analysis.json").read_text())
        assert report["summary"]["total"] > 0

    def test_hdf5_contains_expected_datasets(self, su2_case, tmp_path):
        h5py = pytest.importorskip("h5py")
        out = tmp_path / "case.cfdx.h5"
        convert(su2_case, output=out)
        with h5py.File(out, "r") as f:
            assert f.attrs["source_solver"] == "SU2"
            assert "points" in f
            assert "face_vertices" in f
            assert f.attrs["mesh_topology"] == "cfdx-csr-v1"
            assert "fields" in f
            assert "case_setup_json" in f.attrs

    def test_cfdx_file_is_normalized(self, su2_case, tmp_path):
        first = tmp_path / "case.cfdx.h5"
        convert(su2_case, output=first)

        second = tmp_path / "normalized.cfdx.h5"
        result = convert(first, output=second)

        assert second.exists()
        assert result.success is True
        assert result.has_mesh is True
        assert len(result.fields) > 0
        # provenance of the original solver format is preserved
        assert result.source.solver == "SU2"
        assert result.source.format == "su2"

    def test_get_converter_cfdx_reads_back(self, su2_case, tmp_path):
        out = tmp_path / "case.cfdx.h5"
        original = convert(su2_case, output=out)
        reread = get_converter("cfdx").convert(out)

        assert reread.field_names == original.field_names
        assert reread.has_mesh is True
        assert np.asarray(reread.mesh["points"]).shape[1] == 3

    def test_cfdx_writer_without_mesh_still_writes_setup(self, tmp_path):
        result = ConversionResult(
            source=SourceInfo(solver="Fluent"),
            case=CaseSetup(),
        )
        out = tmp_path / "empty.cfdx.h5"
        from cfdx.io.converter import write_result

        write_result(result, out)
        assert out.exists()


# --------------------------------------------------------------------------
# Gap analysis: nothing is silently dropped
# --------------------------------------------------------------------------


class TestGapAnalysis:
    def test_missing_mesh_is_recorded(self, su2_case):
        result = convert(su2_case)
        assert result.has_mesh is True
        assert not any(
            f.category.startswith("mesh") and f.severity == Severity.UNAVAILABLE
            for f in result.gap_analysis.findings
        )

    def test_missing_materials_are_recorded(self, su2_case):
        result = convert(su2_case)
        assert result.materials == []
        assert any(
            f.category == "material" and f.severity == Severity.UNAVAILABLE
            for f in result.gap_analysis.findings
        )

    def test_missing_fields_are_recorded(self, openfoam_case):
        result = convert(openfoam_case)
        assert result.fields == []
        assert any(
            f.severity == Severity.UNAVAILABLE
            and "field" in f.category + f.feature
            for f in result.gap_analysis.findings
        )

    def test_unreadable_openfoam_field_is_reported(self, openfoam_case):
        zero = openfoam_case / "0"
        zero.mkdir()
        (zero / "k").write_text(
            _foam_header("k", "volScalarField")
            + "internalField nonuniform List<volScalar> 2 ((0 0))\n;\n"
        )
        result = convert(openfoam_case)
        assert "k" not in [name for name, _ in result.fields]
        assert any(
            f.feature == "k" and f.severity == Severity.UNSUPPORTED_NONBLOCK
            for f in result.gap_analysis.findings
        )

    def test_missing_config_is_reported(self, tmp_path):
        mesh = tmp_path / "case.su2"
        shutil.copy(
            os.path.join(DATA_DIR, "su2", "mesh_NACA0012_inv.su2"), mesh
        )
        result = convert(mesh)
        assert result.has_mesh is True
        assert any(
            f.category == "config" for f in result.gap_analysis.findings
        )

    def test_openfoam_without_polymesh_is_blocking(self, tmp_path):
        case = tmp_path / "case"
        case.mkdir()
        result = convert(case)
        assert result.success is False
        assert result.gap_analysis.n_unsupported_blocking() > 0
        assert not result.has_mesh

    def test_result_without_any_data_records_all_gaps(self):
        result = ConversionResult()
        from cfdx.io.converter import _record_missing_data

        _record_missing_data(result)
        categories = {f.category for f in result.gap_analysis.findings}
        assert {"mesh", "fields", "material", "boundary"} <= categories
        assert result.success is True  # missing data is not blocking

    def test_findings_are_not_duplicated(self, su2_case):
        result = convert(su2_case)
        findings = result.gap_analysis.findings
        assert len(findings) == len({(f.category, f.feature) for f in findings})
