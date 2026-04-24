import wickd as w


def setup():
    w.reset_space()
    w.add_space("o", "fermion", "occupied", ["i", "j"])
    w.add_space("a", "fermion", "general", ["u", "v"])
    w.add_space("v", "fermion", "unoccupied", ["a", "b", "c"])


def test_tensor():
    """Test the Tensor class"""
    setup()

    # Create a tensor from a string
    t = w.tensor("T^{v0}_{a0}", w.sym.none)
    assert str(t) == "T^{v0}_{a0}[N]"

    # Create a tensor from a string
    t = w.tensor("T^{o_0,o_1}_{v_0,v_1}", w.sym.none)
    assert str(t) == "T^{o0,o1}_{v0,v1}[N]"

    # Create a tensor with one lower index
    t = w.tensor("T", ["o_0"], [], w.sym.none)
    assert str(t) == "T^{}_{o0}[N]"

    # Create a tensor with one upper index
    t = w.tensor("T", [], ["o_0"], w.sym.none)
    assert str(t) == "T^{o0}_{}[N]"

    # Create a tensor with lower and upper indices
    t = w.tensor("T", ["v_0"], ["o_0"], w.sym.none)
    assert str(t) == "T^{o0}_{v0}[N]"

    # Create a tensor with lower and upper indices
    t = w.tensor("T", ["v_0", "a_0"], ["o_0", "a_1"], w.sym.none)
    assert str(t) == "T^{o0,a1}_{v0,a0}[N]"
    assert t.latex() == "{T}^{i v}_{a u}"


def test_tensor_symmetry_types():
    """Test all three symmetry types via constructor"""
    setup()

    t_anti = w.tensor("T", ["v0"], ["o0"], w.sym.anti)
    t_symm = w.tensor("T", ["v0"], ["o0"], w.sym.symm)
    t_none = w.tensor("T", ["v0"], ["o0"], w.sym.none)

    assert t_anti.symmetry() == w.sym.anti
    assert t_symm.symmetry() == w.sym.symm
    assert t_none.symmetry() == w.sym.none

    assert not t_anti.is_complex_conjugate()
    assert not t_symm.is_complex_conjugate()
    assert not t_none.is_complex_conjugate()

    # str() output: anti emits no modifier, others do
    assert str(t_anti) == "T^{o0}_{v0}"
    assert str(t_symm) == "T^{o0}_{v0}[S]"
    assert str(t_none) == "T^{o0}_{v0}[N]"


def test_tensor_modifier_syntax_all_six():
    """Test all 6 symmetry × conjugation combinations via string modifier syntax"""
    setup()

    cases = [
        # (input_str,          expected_sym,   expected_conj, expected_str)
        ("T^{v0}_{o0}",        w.sym.anti,     False,         "T^{v0}_{o0}"),
        ("T^{v0}_{o0}[A]",     w.sym.anti,     False,         "T^{v0}_{o0}"),
        ("T^{v0}_{o0}[A,*]",   w.sym.anti,     True,          "T^{v0}_{o0}[A,*]"),
        ("T^{v0}_{o0}[S]",     w.sym.symm,     False,         "T^{v0}_{o0}[S]"),
        ("T^{v0}_{o0}[S,*]",   w.sym.symm,     True,          "T^{v0}_{o0}[S,*]"),
        ("T^{v0}_{o0}[N]",     w.sym.none,     False,         "T^{v0}_{o0}[N]"),
        ("T^{v0}_{o0}[N,*]",   w.sym.none,     True,          "T^{v0}_{o0}[N,*]"),
    ]

    for s, expected_sym, expected_conj, expected_str in cases:
        t = w.tensor(s, w.sym.anti)  # default symmetry should be overridden by modifier
        assert t.symmetry() == expected_sym, f"symmetry mismatch for {s!r}"
        assert t.is_complex_conjugate() == expected_conj, f"conjugate mismatch for {s!r}"
        assert str(t) == expected_str, f"str() mismatch for {s!r}"


def test_tensor_modifier_overrides_default_symmetry():
    """Modifier in string overrides the symmetry= argument"""
    setup()

    # Pass sym.none as default but modifier says [A] → should be anti
    t = w.tensor("T^{v0}_{o0}[A]", w.sym.none)
    assert t.symmetry() == w.sym.anti

    # Pass sym.anti as default but modifier says [S] → should be symm
    t = w.tensor("T^{v0}_{o0}[S]", w.sym.anti)
    assert t.symmetry() == w.sym.symm

    # Pass sym.symm as default but modifier says [N] → should be none
    t = w.tensor("T^{v0}_{o0}[N]", w.sym.symm)
    assert t.symmetry() == w.sym.none


def test_tensor_equality_symmetry_distinguishes():
    """Tensors with different symmetry types are not equal"""
    setup()

    t_anti = w.tensor("T^{v0}_{o0}", w.sym.anti)
    t_symm = w.tensor("T^{v0}_{o0}[S]", w.sym.anti)
    t_none = w.tensor("T^{v0}_{o0}[N]", w.sym.anti)

    assert t_anti == t_anti
    assert t_symm == t_symm
    assert t_none == t_none

    assert not (t_anti == t_symm)
    assert not (t_anti == t_none)
    assert not (t_symm == t_none)


def test_tensor_equality_conjugation_distinguishes():
    """Tensors differing only in conjugation are not equal"""
    setup()

    t     = w.tensor("T^{v0}_{o0}",      w.sym.anti)
    t_cc  = w.tensor("T^{v0}_{o0}[A,*]", w.sym.anti)
    t_s   = w.tensor("T^{v0}_{o0}[S]",   w.sym.anti)
    t_s_cc = w.tensor("T^{v0}_{o0}[S,*]", w.sym.anti)

    assert t != t_cc
    assert t_s != t_s_cc
    assert t != t_s_cc


def test_tensor_str_round_trip():
    """str(t) round-trips back to an equal tensor when parsed"""
    setup()

    originals = [
        w.tensor("T^{v0}_{o0}",       w.sym.anti),
        w.tensor("T^{v0}_{o0}[A,*]",  w.sym.anti),
        w.tensor("T^{v0}_{o0}[S]",    w.sym.anti),
        w.tensor("T^{v0}_{o0}[S,*]",  w.sym.anti),
        w.tensor("T^{v0}_{o0}[N]",    w.sym.anti),
        w.tensor("T^{v0}_{o0}[N,*]",  w.sym.anti),
    ]

    for t in originals:
        t2 = w.tensor(str(t), w.sym.anti)
        assert t == t2, f"round-trip failed for {str(t)!r}"


def test_tensor_property_accessors():
    """label(), lower(), upper() return correct values"""
    setup()

    t = w.tensor("T", ["v0", "o0"], ["a0"], w.sym.symm)
    assert t.label() == "T"
    assert len(t.lower()) == 2
    assert len(t.upper()) == 1
    assert t.symmetry() == w.sym.symm
    assert not t.is_complex_conjugate()

    t_cc = w.tensor("T^{a0}_{v0,o0}[S,*]", w.sym.anti)
    assert t_cc.label() == "T"
    assert t_cc.symmetry() == w.sym.symm
    assert t_cc.is_complex_conjugate()


def test_tensor_label_and_indices():
    """Verify label and index accessors for each symmetry type"""
    setup()

    for sym_val, modifier in [(w.sym.anti, ""), (w.sym.symm, "[S]"), (w.sym.none, "[N]")]:
        t = w.tensor("T^{v0}_{o0}" + modifier, w.sym.anti)
        assert t.label() == "T"
        assert t.symmetry() == sym_val


def test_tensor_lt():
    """__lt__ imposes a total order consistent with __eq__"""
    setup()

    t_a = w.tensor("A^{v0}_{o0}", w.sym.anti)
    t_b = w.tensor("B^{v0}_{o0}", w.sym.anti)
    assert t_a < t_b
    assert not (t_b < t_a)
    assert not (t_a < t_a)


def test_tensor_rank_and_indices():
    """rank() and indices() are consistent with lower/upper"""
    setup()

    t = w.tensor("T", ["v0", "o0"], ["a0"], w.sym.anti)
    assert t.rank() == 3  # 2 lower + 1 upper
    idxs = t.indices()
    assert len(idxs) == 3  # all distinct

    # scalar tensor
    t0 = w.tensor("E", [], [], w.sym.none)
    assert t0.rank() == 0
    assert t0.indices() == []


def test_tensor_signature_and_symmetry_factor():
    """signature() and symmetry_factor() return sensible values"""
    setup()

    # two lower in occupied, one upper in virtual → signature entry for o-space
    t = w.tensor("T", ["o0", "o1"], ["v0"], w.sym.anti)
    sig = t.signature()
    assert isinstance(sig, list)
    assert len(sig) == w.num_spaces()  # one entry per space

    sf = t.symmetry_factor()
    # two occupied lower indices → 2!  = 2; one upper → 1! = 1 → factor = 2
    assert sf == 2


def test_tensor_adjoint():
    """adjoint() toggles is_complex_conjugate, leaves label/indices intact"""
    setup()

    t = w.tensor("T^{v0}_{o0}", w.sym.anti)
    assert not t.is_complex_conjugate()
    t_adj = t.adjoint()
    assert t_adj.is_complex_conjugate()
    assert t_adj.label() == t.label()
    assert t_adj.lower() == t.lower()
    assert t_adj.upper() == t.upper()
    # adjoint of adjoint is the original
    assert t_adj.adjoint() == t
