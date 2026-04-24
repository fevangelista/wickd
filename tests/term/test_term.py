import wickd as w


def test_term():
    """Test the SymbolicTerm and Term classes"""
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])
    w.add_space("a", "fermion", "general", ["u", "v"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])

    # Create an empty term
    term = w.Term()
    assert str(term) == "1"

    # Create a term from a symbolic term
    term = w.SymbolicTerm()
    term.set_normal_ordered(True)
    term.add([w.cre("v_0"), w.ann("o_0")])
    term = w.Term(term)
    assert str(term) == "{ a+(v0) a-(o0) }"

    # Create a term from scratch
    term = w.Term()
    term.set_normal_ordered(True)
    term.set(w.rational(1, 2))
    term.add([w.cre("v_0"), w.ann("o_0")])
    term.add(w.tensor("T", ["v_0"], ["o_0"], w.sym.anti))
    assert str(term) == "1/2 T^{o0}_{v0} { a+(v0) a-(o0) }"
    assert term.latex() == "+\\frac{1}{2} {T}^{i}_{a} \\{ \\hat{a}^{a}_{i} \\}"


def _setup():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b"])


def test_symbolicterm_eq_lt():
    """__eq__ and __lt__ on SymbolicTerm"""
    _setup()

    t1 = w.SymbolicTerm()
    t1.set_normal_ordered(True)
    t1.add(w.cre("v_0"))

    t2 = w.SymbolicTerm()
    t2.set_normal_ordered(True)
    t2.add(w.cre("v_0"))

    t3 = w.SymbolicTerm()
    t3.set_normal_ordered(True)
    t3.add(w.cre("v_1"))

    assert t1 == t2
    assert not (t1 == t3)
    # __lt__ gives a total order; exactly one of (t1 < t3) or (t3 < t1) is True
    assert (t1 < t3) != (t3 < t1)
    assert not (t1 < t1)


def test_symbolicterm_nops():
    """nops() counts the SQ operators"""
    _setup()

    t = w.SymbolicTerm()
    assert t.nops() == 0
    t.add(w.cre("v_0"))
    assert t.nops() == 1
    t.add(w.ann("o_0"))
    assert t.nops() == 2


def test_symbolicterm_normal_order_predicates():
    """is_vacuum_normal_ordered, is_labeled_normal_ordered, is_creation_then_annihilation"""
    _setup()

    # a+(v0) a-(o0): cre before ann → vacuum + labeled normal ordered
    t_no = w.SymbolicTerm()
    t_no.set_normal_ordered(True)
    t_no.add([w.cre("v_0"), w.ann("o_0")])
    assert t_no.is_vacuum_normal_ordered()
    assert t_no.is_labeled_normal_ordered()
    assert t_no.is_creation_then_annihilation()

    # a-(o0) a+(v0): ann before cre → NOT labeled normal ordered
    t_rev = w.SymbolicTerm()
    t_rev.set_normal_ordered(False)
    t_rev.add([w.ann("o_0"), w.cre("v_0")])
    assert not t_rev.is_labeled_normal_ordered()
    assert not t_rev.is_creation_then_annihilation()


def test_symbolicterm_adjoint():
    """adjoint() swaps cre↔ann on every operator"""
    _setup()

    t = w.SymbolicTerm()
    t.set_normal_ordered(True)
    t.add([w.cre("v_0"), w.ann("o_0")])
    t_adj = t.adjoint()
    ops = t_adj.ops()
    # Adjoint reverses order and swaps type
    assert ops[0].type() == w.type.cre  # ann("o_0") → cre("o_0")
    assert ops[1].type() == w.type.ann  # cre("v_0") → ann("v_0")


def test_symbolicterm_vacuum_normal_order():
    """vacuum_normal_order() returns a list of (SymbolicTerm, rational) pairs"""
    _setup()

    # Already normal-ordered term should map to itself with coeff 1
    t = w.SymbolicTerm()
    t.set_normal_ordered(False)
    t.add([w.cre("v_0"), w.ann("o_0")])
    result = t.vacuum_normal_order()
    assert isinstance(result, list)
    assert len(result) >= 1
    for sterm, coeff in result:
        assert isinstance(sterm, w.SymbolicTerm)


def test_symbolicterm_compile():
    """compile() returns a non-empty string"""
    _setup()

    t = w.SymbolicTerm()
    t.set_normal_ordered(True)
    t.add(w.tensor("T", ["v_0"], ["o_0"], w.sym.anti))
    t.add([w.cre("v_0"), w.ann("o_0")])
    out = t.compile("einsum")
    assert isinstance(out, str)


def test_term_coefficient_and_symterm():
    """coefficient() and symterm() are accessible and consistent"""
    _setup()

    t = w.Term()
    t.set_normal_ordered(True)
    t.set(w.rational(1, 3))
    t.add([w.cre("v_0"), w.ann("o_0")])

    assert t.coefficient() == w.rational(1, 3)
    st = t.symterm()
    assert isinstance(st, w.SymbolicTerm)
    assert st.nops() == 2
