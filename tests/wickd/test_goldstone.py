import pytest
import wickd as w


@pytest.fixture(autouse=True)
def ov_spaces():
    """Occupied + virtual spaces used by all tests in this module."""
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j", "k", "l"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c", "d"])
    yield
    w.reset_space()


def _contract(ops):
    wt = w.WickTheorem()
    return wt.contract(ops, minrank=0, maxrank=0)


# ── space_type binding ────────────────────────────────────────────────────────


def test_space_type():
    osi = w.osi()
    assert osi.space_type(0) == "occupied"
    assert osi.space_type(1) == "unoccupied"


# ── expression_to_tikz ────────────────────────────────────────────────────────


def test_returns_string():
    F = w.utils.gen_op("f", 1, "ov", "ov")
    T1 = w.op("t", ["v+ o"])
    expr = _contract(F @ T1)
    tikz = w.expression_to_tikz(expr)
    assert isinstance(tikz, str)
    assert len(tikz) > 0


def test_tikzpicture_tags():
    F = w.utils.gen_op("f", 1, "ov", "ov")
    T1 = w.op("t", ["v+ o"])
    expr = _contract(F @ T1)
    tikz = w.expression_to_tikz(expr)
    n_terms = len(list(expr))
    assert tikz.count(r"\begin{tikzpicture}") == n_terms
    assert tikz.count(r"\end{tikzpicture}") == n_terms


def test_standalone_wraps_document():
    F = w.utils.gen_op("f", 1, "ov", "ov")
    T1 = w.op("t", ["v+ o"])
    expr = _contract(F @ T1)
    tikz = w.expression_to_tikz(expr, standalone=True)
    assert r"\documentclass" in tikz
    assert r"\begin{document}" in tikz
    assert r"\end{document}" in tikz
    assert r"\usetikzlibrary{arrows.meta}" in tikz


def test_operator_bars_present():
    F = w.utils.gen_op("f", 1, "ov", "ov")
    T1 = w.op("t", ["v+ o"])
    expr = _contract(F @ T1)
    tikz = w.expression_to_tikz(expr)
    # Each diagram should have draw commands for the operator bars
    assert r"\draw[thick]" in tikz


def test_arrows_present():
    F = w.utils.gen_op("f", 1, "ov", "ov")
    T1 = w.op("t", ["v+ o"])
    expr = _contract(F @ T1)
    tikz = w.expression_to_tikz(expr)
    assert "Stealth" in tikz


def test_two_body_interaction():
    """Two-body Hamiltonian contracted with T2 — multiple lines per diagram."""
    V = w.op("v", ["o+ o+ v v"])
    T2 = w.op("t", ["v+ v+ o o"])
    expr = _contract(V @ T2)
    tikz = w.expression_to_tikz(expr)
    assert isinstance(tikz, str)
    assert r"\begin{tikzpicture}" in tikz


def test_empty_expression():
    """An expression with no fully contracted terms yields an empty string."""
    expr = w.Expression()
    tikz = w.expression_to_tikz(expr)
    assert tikz == ""
