"""Reference methods stay distinct from native production availability."""

from cfdx.modern_methods import (
    METHODS,
    method_by_key,
    method_families,
    methods_for_family,
)


def test_catalog_contains_traceable_nonlinear_research_candidates():
    assert len({item.key for item in METHODS}) == len(METHODS)
    for key in ("NEWTON_KRYLOV", "ANDERSON", "NGMRES"):
        method = method_by_key(key)
        assert (
            "no native" in method.integration.casefold()
            and "certified" in method.integration.casefold()
        )
        assert method.reference.startswith("https://")
        assert method.qualification


def test_adapter_mapping_does_not_imply_numerical_certification():
    for key in ("SIMPLE", "SIMPLEC", "COUPLED", "PCD", "AMG", "WLS", "TVD", "FGMRES"):
        method = method_by_key(key)
        assert method.integration == "Adapter mapping exists"
        assert method.purpose
    assert "Nonlinear acceleration" in method_families()
    assert len(methods_for_family("Nonlinear acceleration")) == 3


def test_search_and_compare_are_read_only_and_explicit():
    import pytest
    from cfdx.modern_methods import compare_methods, search_methods

    assert {m.key for m in search_methods("krylov")} == {"NEWTON_KRYLOV", "FGMRES"}
    assert {m.key for m in search_methods("schur", integration="mapped")} == {
        "COUPLED",
        "PCD",
    }
    assert {m.key for m in search_methods(integration="research")} == {
        "NEWTON_KRYLOV",
        "ANDERSON",
        "NGMRES",
    }
    assert search_methods("xxxxx") == ()
    with pytest.raises(ValueError, match="invalid method integration"):
        search_methods(integration="unknown")
    with pytest.raises(ValueError, match="different numerical methods"):
        compare_methods("AMG", "AMG")
    comparison = compare_methods("COUPLED", "NEWTON_KRYLOV")
    assert "READ ONLY" in comparison
    assert "not certification" in comparison
    assert "no native integration certified" in comparison
    assert "identical case" in comparison
