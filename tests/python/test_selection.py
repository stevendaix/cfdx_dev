"""Headless tests for the shared selection identity scheme."""

import pytest


def test_patch_selection_round_trip() -> None:
    from cfdx.selection import parse_selection_id, patch_selection_id

    stable_id = patch_selection_id("wall")
    assert stable_id == "patch:wall"
    target = parse_selection_id(stable_id)
    assert target.kind == "patch"
    assert target.stable_id == stable_id
    assert target.name == "wall"
    assert target.region is None
    assert target.face_index is None


def test_region_scoped_patch_selection_round_trip() -> None:
    from cfdx.selection import parse_selection_id, patch_selection_id

    stable_id = patch_selection_id("inlet", region="fluid")
    assert stable_id == "region:fluid/patch:inlet"
    target = parse_selection_id(stable_id)
    assert target.kind == "patch"
    assert target.name == "inlet"
    assert target.region == "fluid"


def test_region_selection_round_trip() -> None:
    from cfdx.selection import parse_selection_id, region_selection_id

    stable_id = region_selection_id("fluid")
    assert stable_id == "region:fluid"
    target = parse_selection_id(stable_id)
    assert target.kind == "region"
    assert target.name == "fluid"


def test_face_selection_round_trip() -> None:
    from cfdx.selection import face_selection_id, parse_selection_id

    stable_id = face_selection_id("region:fluid/patch:inlet", 3)
    assert stable_id == "region:fluid/patch:inlet/face:3"
    target = parse_selection_id(stable_id)
    assert target.kind == "face"
    assert target.name == "inlet"
    assert target.region == "fluid"
    assert target.face_index == 3


def test_default_region_patch_has_no_prefix() -> None:
    from cfdx.selection import patch_selection_id

    assert patch_selection_id("wall", region="default") == "patch:wall"
    assert patch_selection_id("wall", region=None) == "patch:wall"


def test_mesh_model_contract_ids_match_the_selection_scheme() -> None:
    from cfdx.mesh_model import MeshPatch, MeshRegion
    from cfdx.selection import patch_selection_id, region_selection_id

    patch = MeshPatch(0, "inlet", "inlet", 3, "fluid")
    assert patch.stable_id == patch_selection_id("inlet", region="fluid")
    plain = MeshPatch(1, "wall", "wall", 4)
    assert plain.stable_id == patch_selection_id("wall")
    region = MeshRegion("fluid", 8, 19, 3, (patch,), ("c0",))
    assert region.stable_id == region_selection_id("fluid")


def test_selection_labels_are_human_readable() -> None:
    from cfdx.selection import parse_selection_id, selection_label

    def label(stable_id):
        return selection_label(parse_selection_id(stable_id))

    assert label("region:fluid") == "Region fluid"
    assert label("patch:wall") == "Patch wall"
    assert label("region:fluid/patch:wall") == "Patch wall — region fluid"
    assert label("patch:inlet/face:3") == "Face 3 on patch inlet"
    scoped_face = label("region:fluid/patch:inlet/face:3")
    assert scoped_face == "Face 3 on patch inlet — region fluid"


@pytest.mark.parametrize(
    "stable_id",
    [
        "",
        " patch:inlet",
        "patch:inlet ",
        "patch:",
        "region:",
        "cell:12",
        "patch:inlet/face:03",
        "patch:inlet/face:-1",
        "patch:inlet/face:1/face:2",
        "region:fluid/region:solid",
        "region:fluid/patch:wall/face:",
        "patch:inlet/patch:wall",
        "region:fluid/",
        "face:3",
    ],
)
def test_parse_selection_id_rejects_invalid_ids(stable_id) -> None:
    from cfdx.selection import parse_selection_id

    with pytest.raises(ValueError, match="invalid selection stable id"):
        parse_selection_id(stable_id)


def test_face_selection_id_rejects_invalid_inputs() -> None:
    from cfdx.selection import face_selection_id

    with pytest.raises(ValueError, match="invalid selection stable id"):
        face_selection_id("patch:inlet", -1)
    with pytest.raises(ValueError, match="invalid selection stable id"):
        face_selection_id("patch:inlet", True)
    with pytest.raises(ValueError, match="invalid selection stable id"):
        face_selection_id("region:fluid", 0)
    with pytest.raises(ValueError, match="invalid selection stable id"):
        face_selection_id("patch:inlet/face:1", 0)


def test_builders_reject_non_canonical_names() -> None:
    from cfdx.selection import patch_selection_id, region_selection_id

    with pytest.raises(ValueError, match="invalid selection stable id"):
        patch_selection_id(" wall")
    with pytest.raises(ValueError, match="invalid selection stable id"):
        patch_selection_id("wall", region="fluid/solid")
    with pytest.raises(ValueError, match="invalid selection stable id"):
        region_selection_id("")
