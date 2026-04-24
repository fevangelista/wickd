import wickd as w


def test_sqoperator():
    """Test the SQOperator class"""
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])
    w.add_space("p", "boson", "unoccupied", ["u", "v"])

    cop = w.sqoperator("o_0", w.type.cre)
    assert str(cop) == "a+(o0)"
    cop = w.sqoperator("v_1", w.type.ann)
    assert str(cop) == "a-(v1)"
    cop = w.sqoperator("p_0", w.type.cre)
    assert str(cop) == "b+(p0)"
    cop = w.sqoperator("p_1", w.type.ann)
    assert str(cop) == "b-(p1)"

    cop = w.cre("o_0")
    assert str(cop) == "a+(o0)"
    aop = w.ann("v_1")
    assert str(aop) == "a-(v1)"


def test_sqoperator2():
    """Test the SQOperator class"""
    w.reset_space()
    w.add_space("c", "fermion", "occupied", ["i", "j"])
    w.add_space("a", "fermion", "general", ["u", "v"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])
    assert w.cre("c_0") < w.ann("c_0")
    assert w.cre("a_0") < w.ann("c_0")
    assert w.cre("v_0") < w.ann("c_0")
    assert w.cre("v_0") < w.ann("a_0")
    assert w.cre("c_0") < w.cre("c_1")
    assert w.cre("c_0") < w.cre("a_1")
    assert w.ann("c_1") < w.ann("c_0")


def test_sqoperator3():
    """Test the SQOperator class"""
    w.reset_space()
    w.add_space("c", "fermion", "occupied", ["i", "j"])
    w.add_space("a", "fermion", "general", ["u", "v"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])

    assert w.ann("c_0").normal_ordered_less(w.cre("c_0"))
    assert not w.cre("c_0").normal_ordered_less(w.cre("c_0"))
    assert w.cre("c_1").normal_ordered_less(w.cre("c_0"))
    assert not w.cre("c_0").normal_ordered_less(w.cre("c_1"))
    assert not w.ann("c_0").normal_ordered_less(w.ann("c_0"))
    assert w.ann("c_0").normal_ordered_less(w.ann("c_1"))
    assert not w.ann("c_1").normal_ordered_less(w.ann("c_0"))

    assert w.cre("v_0").normal_ordered_less(w.ann("v_0"))
    assert not w.cre("v_0").normal_ordered_less(w.cre("v_0"))
    assert w.cre("v_0").normal_ordered_less(w.cre("v_1"))
    assert not w.cre("v_1").normal_ordered_less(w.cre("v_0"))
    assert not w.ann("v_0").normal_ordered_less(w.ann("v_0"))
    assert w.ann("v_1").normal_ordered_less(w.ann("v_0"))
    assert not w.ann("v_0").normal_ordered_less(w.ann("v_1"))

    assert w.cre("a_0").normal_ordered_less(w.ann("a_0"))
    assert not w.cre("a_0").normal_ordered_less(w.cre("a_0"))
    assert w.cre("a_0").normal_ordered_less(w.cre("a_1"))
    assert not w.cre("a_1").normal_ordered_less(w.cre("a_0"))
    assert not w.ann("a_0").normal_ordered_less(w.ann("a_0"))
    assert w.ann("a_1").normal_ordered_less(w.ann("a_0"))
    assert not w.ann("a_0").normal_ordered_less(w.ann("a_1"))

    assert w.ann("c_0").normal_ordered_less(w.ann("v_0"))
    assert not w.ann("v_0").normal_ordered_less(w.ann("c_0"))


def test_sqoperator_new_accessors():
    """Test accessors added to complete the Python API"""
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])   # space index 0
    w.add_space("v", "fermion", "unoccupied", ["a", "b"]) # space index 1

    cop = w.cre("v_0")
    aop = w.ann("o_0")

    # space() returns the integer space index
    assert cop.space() == 1  # v is space 1
    assert aop.space() == 0  # o is space 0

    # is_creation()
    assert cop.is_creation()
    assert not aop.is_creation()

    # is_quasiparticle_creation(): cre of unoccupied = particle creation = True
    #                               ann of occupied  = hole creation      = True
    assert cop.is_quasiparticle_creation()
    assert aop.is_quasiparticle_creation()
    assert not w.ann("v_0").is_quasiparticle_creation()
    assert not w.cre("o_0").is_quasiparticle_creation()

    # op_symbol() returns the field-type base symbol: "a" for fermions, "b" for bosons
    assert cop.op_symbol() == "a"
    assert aop.op_symbol() == "a"

    # adjoint(): cre → ann, ann → cre; index preserved
    adj = cop.adjoint()
    assert adj.type() == w.type.ann
    assert adj.index() == cop.index()
    adj2 = aop.adjoint()
    assert adj2.type() == w.type.cre
    assert adj2.index() == aop.index()

    # compile() is exposed but not implemented (raises RuntimeError)
    import pytest
    with pytest.raises(RuntimeError):
        cop.compile("einsum")
